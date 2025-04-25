import unittest
import glob
import numpy as np

from stcrpy.tcr_processing import TCRParser
from stcrpy.tcr_geometry import TCRDock, TCRCoM, TCRGeom


class TestTCRGeometry(unittest.TestCase):
    def test_TCRDock_init(self):
        parser = TCRParser.TCRParser()
        pdb_file = "./test_files/5hyj.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)
        [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]

    def test_calculate_docking_angle_5hyj(self):
        parser = TCRParser.TCRParser()
        pdb_file = "./test_files/5hyj.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)

        tcr_docks = [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]
        assert all(
            [
                x.calculate_docking_angle() < 50.0
                and x.calculate_docking_angle() > 40.0
                for x in tcr_docks
            ]
        )

    def test_calculate_docking_angle_7l1d(self):
        parser = TCRParser.TCRParser()
        pdb_file = "./test_files/7l1d.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)

        tcr_docks = [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]
        assert all(
            [
                x.calculate_docking_angle() < 50.0
                and x.calculate_docking_angle() > 40.0
                for x in tcr_docks
            ]
        )

    def test_calculate_docking_angle_7rrg(self):
        parser = TCRParser.TCRParser()
        pdb_file = "./test_files/7rrg.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)

        tcr_docks = [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]
        assert all(
            [
                x.calculate_docking_angle() < 80.0
                and x.calculate_docking_angle() > 70.0
                for x in tcr_docks
            ]
        )

    def test_calculate_docking_angle_of_docks(self):
        parser = TCRParser.TCRParser()
        dock_pdb_files = glob.glob(
            "./test_files/TCRHaddock_test_files/387937-tcr_6eqa_mel5_bulged/structures/it1/renumbered_complex*.pdb"
        )
        dock_pdb_files.sort()
        for pdb_file in dock_pdb_files:
            tcr = parser.get_tcr_structure("test", pdb_file)
            [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]

    def test_docking_angle_reverse_docks(self):
        parser = TCRParser.TCRParser()
        reverse_dock_108 = "./test_files/aligned_complex_108.pdb"
        dock_63 = "./test_files/aligned_complex_63.pdb"
        dock_pdb_files = [dock_63, reverse_dock_108]
        docking_angles = []
        for pdb_file in dock_pdb_files:
            tcr = parser.get_tcr_structure("test", pdb_file)

            tcr_docks = [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]
            docking_angles.append(tcr_docks[0].calculate_docking_angle())

        print(docking_angles)

    def test_MH1_TCRCoM(self):
        parser = TCRParser.TCRParser()
        pdb_file = "./test_files/4nhu.pdb"
        tcr_structure = parser.get_tcr_structure("test", pdb_file)
        tcr_com = TCRCoM.MHCI_TCRCoM()
        # TEST DEPRECATED
        # for tcr in tcr_structure.get_TCRs():
        #     r, theta, phi = tcr_com.calculate_geometry(
        #         tcr,
        #         save_aligned_as=f"./test_files/out/aligned_test_{tcr.id}.pdb",
        #     )
        #     print(r, theta, phi)
        # pdb_file = "./stcrpy/stcrpy/tcr_geometry/reference_data/dock_reference_1_imgt_numbered.pdb"
        # tcr_structure = parser.get_tcr_structure("test", pdb_file)
        # tcr_com = TCRCoM.MHCI_TCRCoM()
        # for tcr in tcr_structure.get_TCRs():
        #     r, theta, phi = tcr_com.calculate_geometry(
        #         tcr,
        #         save_aligned_as=f"./test_files/out/aligned_test_dock_ref_mhcI_{tcr.id}.pdb",
        #     )
        #     print(r, theta, phi)

    def test_MH2_TCRCoM(self):
        parser = TCRParser.TCRParser()
        pdb_file = "./test_files/6r0e.cif"
        tcr_structure = parser.get_tcr_structure("test", pdb_file)
        tcr_com = TCRCoM.MHCII_TCRCoM()
        # TEST DEPRECATED
        # for tcr in tcr_structure.get_TCRs():
        #     r, theta, phi = tcr_com.calculate_geometry(
        #         tcr,
        #         save_aligned_as=f"./test_files/out/aligned_test_{tcr.id}.pdb",
        #     )

        # pdb_file = (
        #     "./stcrpy/stcrpy/tcr_geometry/include/dock_reference_2_imgt_numbered.pdb"
        # )
        # tcr_structure = parser.get_tcr_structure("test", pdb_file)
        # tcr_com = TCRCoM.MHCII_TCRCoM()
        # for tcr in tcr_structure.get_TCRs():
        #     r, theta, phi = tcr_com.calculate_geometry(
        #         tcr,
        #         save_aligned_as=f"./test_files/out/aligned_test_dock_ref_mhcII_{tcr.id}.pdb",
        #     )
        #     print(r, theta, phi)

    def testTCRGeom(self):
        parser = TCRParser.TCRParser()
        pdb_files = glob.glob("./test_files/TCRCoM_test_files/*.cif")
        # 'TCRpy/test/test_files/TCRCoM_test_files/7sg0.cif')
        for file in pdb_files:
            file_id = file.split("/")[-1].split(".")[0]
            print(file_id)
            tcr = parser.get_tcr_structure(file_id, file)
            for x in tcr.get_TCRs():
                try:
                    x.geometry = TCRGeom.TCRGeom(
                        x,
                        save_aligned_as=f"./test_files/out/{file_id}_aligned.pdb",
                    )
                    print(x.geometry)
                except Exception as e:
                    print(e)

    def test_TCR_geom_methods(self):
        parser = TCRParser.TCRParser()
        test_file = "./test_files/8gvb.cif"
        tcr = list(parser.get_tcr_structure("8gvb", test_file).get_TCRs())[0]
        geometry = tcr.calculate_docking_geometry()
        assert "scanning_angle" in geometry

    def test_calculate_docking_angle_cys_method(self):
        import stcrpy

        pdb_files = [
            "./test_files/5hyj.pdb",
            "./test_files/7l1d.pdb",
            "./test_files/7rrg.pdb",
        ]
        tcrs = stcrpy.load_TCRs(pdb_files)

        true_scanning_angles = [
            42.9581,
            47.4101,
            47.4101,
            73.7909,
        ]  # [5hyj, 7l1d, 7rrg] - repeat value is for TCR repeats in PDB file parsed as independent TCR structures
        true_pitch_angles = [12.3062, 3.46555, 3.46555, 12.0141]

        cys_crossing_angles = []
        com_crossing_angles = []
        for i, tcr in enumerate(tcrs):
            cys_crossing_angles.append(tcr.get_scanning_angle(mode="cys"))
            com_crossing_angles.append(tcr.get_scanning_angle(mode="com"))

        self.assertTrue(
            np.sum(
                np.abs(
                    np.asarray(cys_crossing_angles) - np.asarray(true_scanning_angles)
                )
            )
            / len(cys_crossing_angles)
            < 1.5
        )

        # pitch = tcr.get_pitch_angle()
        # self.assertAlmostEqual(pitch, true_pitch_angles[i])

    # def test_calculate_docking_angle_com_method(self):
    #     import stcrpy

    #     pdb_files = [
    #         "./test_files/5hyj.pdb",
    #         "./test_files/7l1d.pdb",
    #         "./test_files/7rrg.pdb",
    #     ]
    #     tcrs = stcrpy.load_TCRs(pdb_files)

    #     true_scanning_angles = [42.9581, 47.4101, 73.7909]  # [5hyj, 7l1d, 7rrg]
    #     true_pitch_angles = [12.3062, 3.46555, 12.0141]
    #     for i, tcr in enumerate(tcrs):
    #         crossing_angle = tcr.get_scanning_angle(mode="com")
    #         self.assertAlmostEqual(
    #             crossing_angle,
    #             true_scanning_angles[i],
    #         )
    #         pitch = tcr.get_pitch_angle()
    #         self.assertAlmostEqual(pitch, true_pitch_angles[i])

    def test_get_alpha_helices(self):
        import stcrpy

        tcrs = stcrpy.load_TCRs(glob.glob("test_files/TCRCoM_test_files/*.cif"))

        for tcr in tcrs:
            if (
                len(tcr.get_MHC()) == 1
                and hasattr(tcr.get_MHC()[0], "MHC_type")
                and tcr.get_MHC()[0].MHC_type in ["MH1", "MH2"]
            ):
                tcr_geom = TCRGeom.TCRGeom(tcr)
                # since tcr should be aligned to reference MHC mhc vector should be approximately [0, 1, 0]
                self.assertAlmostEqual(
                    np.dot(
                        tcr_geom._get_mhc_helix_vectors(tcr.get_MHC()[0]),
                        np.asarray([0, 1, 0]),
                    ),
                    1,
                    places=1,
                )

    def test_rudolph_scanning_angle(self):
        import stcrpy

        tcrs = stcrpy.load_TCRs(
            glob.glob("test_files/TCRGeom_rudolph_test_files/*.cif")
        )

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
