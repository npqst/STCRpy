import unittest
import glob
import numpy as np

import stcrpy
from stcrpy.tcr_processing import TCRParser
from stcrpy.tcr_geometry import TCRDock, TCRGeom

STRUCTURES = "./test_files/structures"
PARSER_FILES = "./test_files/TCRParser_test_files"
RUDOLPH_FILES = "./test_files/TCRGeom_rudolph_test_files"
HADDOCK_FILES = "./test_files/TCRHaddock_test_files"


class TestTCRGeometry(unittest.TestCase):
    def test_TCRDock_init(self):
        tcr = stcrpy.load_TCRs(f"{PARSER_FILES}/5hyj.pdb")
        [TCRDock.TCRDock(x) for x in tcr]

    def test_calculate_docking_angle_5hyj(self):
        tcr = stcrpy.load_TCRs(f"{PARSER_FILES}/5hyj.pdb")
        tcr_docks = [TCRDock.TCRDock(x) for x in tcr]
        assert all(
            x.calculate_docking_angle() < 50.0 and x.calculate_docking_angle() > 40.0
            for x in tcr_docks
        )

    def test_calculate_docking_angle_7l1d(self):
        tcr = stcrpy.load_TCRs(f"{STRUCTURES}/7l1d.cif")
        tcr_docks = [TCRDock.TCRDock(x) for x in tcr]
        assert all(
            x.calculate_docking_angle() < 50.0 and x.calculate_docking_angle() > 40.0
            for x in tcr_docks
        )

    def test_calculate_docking_angle_7rrg(self):
        tcr = stcrpy.load_TCRs(f"{STRUCTURES}/7rrg.cif")
        tcr_docks = [TCRDock.TCRDock(x) for x in tcr]
        assert all(
            x.calculate_docking_angle() < 80.0 and x.calculate_docking_angle() > 70.0
            for x in tcr_docks
        )

    def test_calculate_docking_angle_of_docks(self):
        parser = TCRParser.TCRParser()
        dock_pdb_files = sorted(
            glob.glob(f"{HADDOCK_FILES}/predictions/renumbered_complex*.pdb")
        )
        for pdb_file in dock_pdb_files:
            tcr = parser.get_tcr_structure("test", pdb_file)
            [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]

    def testTCRGeom(self):
        test_files = [
            f"{PARSER_FILES}/5hyj.pdb",
            f"{PARSER_FILES}/6r0e.pdb",
            f"{STRUCTURES}/8gvb.cif",
        ]
        parser = TCRParser.TCRParser()
        for file in test_files:
            file_id = file.split("/")[-1].split(".")[0]
            tcr = parser.get_tcr_structure(file_id, file)
            for x in tcr.get_TCRs():
                x.geometry = TCRGeom.TCRGeom(
                    x,
                    save_aligned_as=f"./test_files/out/{file_id}_aligned.pdb",
                )

    def test_TCR_geom_methods(self):
        parser = TCRParser.TCRParser()
        tcr = list(parser.get_tcr_structure("8gvb", f"{STRUCTURES}/8gvb.cif").get_TCRs())[0]
        geometry = tcr.calculate_docking_geometry()
        assert "scanning_angle" in geometry

    def test_calculate_docking_angle_cys_method(self):
        pdb_files = [
            f"{PARSER_FILES}/5hyj.pdb",
            f"{STRUCTURES}/7l1d.cif",
            f"{STRUCTURES}/7rrg.cif",
        ]
        tcrs = stcrpy.load_TCRs(pdb_files)

        # Expected scanning angles per structure [5hyj×2, 7l1d, 7rrg]
        true_scanning_angles = [42.9581, 47.4101, 47.4101, 73.7909]

        cys_crossing_angles = [tcr.get_scanning_angle(mode="cys") for tcr in tcrs]

        self.assertTrue(
            np.sum(np.abs(np.asarray(cys_crossing_angles) - np.asarray(true_scanning_angles)))
            / len(cys_crossing_angles)
            < 1.5
        )

    def test_calculate_docking_angle_com_method(self):
        pdb_files = [
            f"{PARSER_FILES}/5hyj.pdb",
            f"{STRUCTURES}/7l1d.cif",
            f"{STRUCTURES}/7rrg.cif",
        ]
        tcrs = stcrpy.load_TCRs(pdb_files)

        true_scanning_angles = [42.9581, 47.4101, 73.7909]
        for i, tcr in enumerate(tcrs):
            crossing_angle = tcr.get_scanning_angle(mode="com")
            self.assertAlmostEqual(crossing_angle, true_scanning_angles[i], places=2)

    def test_get_alpha_helices(self):
        test_files = [
            f"{PARSER_FILES}/5hyj.pdb",
            f"{PARSER_FILES}/6r0e.pdb",
            f"{STRUCTURES}/8gvb.cif",
        ]
        tcrs = stcrpy.load_TCRs(test_files)

        for tcr in tcrs:
            if (
                len(tcr.get_MHC()) == 1
                and hasattr(tcr.get_MHC()[0], "MHC_type")
                and tcr.get_MHC()[0].MHC_type in ["MH1", "MH2"]
            ):
                tcr_geom = TCRGeom.TCRGeom(tcr)
                # after alignment to reference MHC, the helix vector should be ~[0, 1, 0]
                self.assertAlmostEqual(
                    np.dot(
                        tcr_geom._get_mhc_helix_vectors(tcr.get_MHC()[0]),
                        np.asarray([0, 1, 0]),
                    ),
                    1,
                    places=1,
                )

    def test_rudolph_scanning_angle(self):
        tcrs = stcrpy.load_TCRs(glob.glob(f"{RUDOLPH_FILES}/*.cif"))

        rudolph_scanning_angles = {
            "2ckb": 22,
            "1g6r": 23,
            "1mwa": 23,
            "1fo0": 41,
            "1nam": 40,
            "1kj2": 31,
            "1bd2": 48,
        }

        for tcr in tcrs:
            tcr_geom = TCRGeom.TCRGeom(tcr, mode="rudolph")
            assert (
                abs(
                    np.degrees(tcr_geom.scanning_angle)
                    - rudolph_scanning_angles[tcr.parent.parent.id]
                )
                < 3.0
            )
