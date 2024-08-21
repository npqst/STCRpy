import unittest
import glob
import numpy as np

from ..STCRpy.tcr_processing import TCRParser
from ..STCRpy.tcr_geometry import TCRDock, TCRAngle, TCRCoM, TCRGeom


class TestTCRGeometry(unittest.TestCase):
    def test_TCRDock_init(self):
        parser = TCRParser.TCRParser()
        pdb_file = "./STCRpy/test/test_files/5hyj.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)
        [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]

    def test_calculate_docking_angle_5hyj(self):
        parser = TCRParser.TCRParser()
        pdb_file = "./STCRpy/test/test_files/5hyj.pdb"
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
        pdb_file = "./STCRpy/test/test_files/7l1d.pdb"
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
        pdb_file = "./STCRpy/test/test_files/7rrg.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)

        tcr_docks = [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]
        assert all(
            [
                x.calculate_docking_angle() < 80.0
                and x.calculate_docking_angle() > 70.0
                for x in tcr_docks
            ]
        )

    def test_calculate_tcr_angles(self):
        parser = TCRParser.TCRParser()
        pdb_file = "./STCRpy/test/test_files/5hyj.pdb"
        tcr = parser.get_tcr_structure("test", pdb_file)
        tcr_angle = TCRAngle.abTCRAngle()
        ED_tcr_angles = tcr_angle.get_angles(tcr[0]["ED"])
        all_tcr_angles = tcr_angle.get_angles(tcr)
        assert all_tcr_angles["ED"] == ED_tcr_angles

    def test_calculate_docking_angle_of_docks(self):
        parser = TCRParser.TCRParser()
        dock_pdb_files = glob.glob(
            """/home/quast/Collaborations/SamuelsLab_Weizmann/docking/results/
N17.2/310569-N17-2_NRAS_rank_0/structures/it1/renumbered_complex_*.pdb"""
        )
        dock_pdb_files.sort()
        for pdb_file in dock_pdb_files:
            tcr = parser.get_tcr_structure("test", pdb_file)
            [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]

    def test_docking_angle_reverse_docks(self):
        parser = TCRParser.TCRParser()
        reverse_dock_108 = "./STCRpy/test/test_files/aligned_complex_108.pdb"
        dock_63 = "./STCRpy/test/test_files/aligned_complex_63.pdb"
        dock_pdb_files = [dock_63, reverse_dock_108]
        docking_angles = []
        for pdb_file in dock_pdb_files:
            tcr = parser.get_tcr_structure("test", pdb_file)

            tcr_docks = [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]
            docking_angles.append(tcr_docks[0].calculate_docking_angle())

        print(docking_angles)

    def test_MH1_TCRCoM(self):
        parser = TCRParser.TCRParser()
        pdb_file = "./STCRpy/test/test_files/4nhu.pdb"
        tcr_structure = parser.get_tcr_structure("test", pdb_file)
        tcr_com = TCRCoM.MHCI_TCRCoM()
        for tcr in tcr_structure.get_TCRs():
            r, theta, phi = tcr_com.calculate_geometry(
                tcr,
                save_aligned_as=f"./STCRpy/test/test_files/out/aligned_test_{tcr.id}.pdb",
            )
            print(r, theta, phi)
        pdb_file = "./STCRpy/STCRpy/tcr_geometry/reference_data/dock_reference_1_imgt_numbered.pdb"
        tcr_structure = parser.get_tcr_structure("test", pdb_file)
        tcr_com = TCRCoM.MHCI_TCRCoM()
        for tcr in tcr_structure.get_TCRs():
            r, theta, phi = tcr_com.calculate_geometry(
                tcr,
                save_aligned_as=f"./STCRpy/test/test_files/out/aligned_test_dock_ref_mhcI_{tcr.id}.pdb",
            )
            print(r, theta, phi)

    def test_MH2_TCRCoM(self):
        parser = TCRParser.TCRParser()
        pdb_file = "./STCRpy/test/test_files/6r0e.cif"
        tcr_structure = parser.get_tcr_structure("test", pdb_file)
        tcr_com = TCRCoM.MHCII_TCRCoM()
        for tcr in tcr_structure.get_TCRs():
            r, theta, phi = tcr_com.calculate_geometry(
                tcr,
                save_aligned_as=f"./STCRpy/test/test_files/out/aligned_test_{tcr.id}.pdb",
            )

        pdb_file = "./STCRpy/STCRpy/tcr_geometry/reference_data/dock_reference_2_imgt_numbered.pdb"
        tcr_structure = parser.get_tcr_structure("test", pdb_file)
        tcr_com = TCRCoM.MHCII_TCRCoM()
        for tcr in tcr_structure.get_TCRs():
            r, theta, phi = tcr_com.calculate_geometry(
                tcr,
                save_aligned_as=f"./STCRpy/test/test_files/out/aligned_test_dock_ref_mhcII_{tcr.id}.pdb",
            )
            print(r, theta, phi)

    def test_error_prone_MHC_I_TCRCoM_examples(self):
        parser = TCRParser.TCRParser()
        pdb_files = glob.glob("./STCRpy/test/test_files/TCRCoM_test_files/*")
        pdb_files = [
            "/home/quast/Projects/STCRDab/Data/entries/8d5q/structure/imgt/8d5q.pdb"
        ]
        tcr_com = TCRCoM.MHCI_TCRCoM()

        for file in pdb_files:
            pdb_id = file.split("/")[-1].split(".")[0]
            print(pdb_id)
            tcr_structure = parser.get_tcr_structure(pdb_id, file)
            for tcr in tcr_structure.get_TCRs():
                try:
                    r, theta, phi = tcr_com.calculate_geometry(
                        tcr,
                        save_aligned_as=f"./STCRpy/test/test_files/out/aligned_test_pdb_id_{tcr.id}.pdb",
                    )
                except AssertionError as e:
                    assert str(e) == "No MHC associated with TCR"
                    assert pdb_id in [
                        "8dfw",
                        "8d5p",
                    ]  # TCRs without MHC that should raise this error

                assert all(not np.isnan(x) for x in (r, theta, phi))

    def testTCRGeom(self):
        parser = TCRParser.TCRParser()
        pdb_files = glob.glob("./STCRpy/test/test_files/TCRCoM_test_files/*.cif")
        # 'TCRpy/test/test_files/TCRCoM_test_files/7sg0.cif')
        for file in pdb_files:
            file_id = file.split("/")[-1].split(".")[0]
            print(file_id)
            tcr = parser.get_tcr_structure(file_id, file)
            for x in tcr.get_TCRs():
                try:
                    x.geometry = TCRGeom.TCRGeom(
                        x,
                        save_aligned_as=f"./STCRpy/test/test_files/out/{file_id}_aligned.pdb",
                    )
                    print(x.geometry)
                except Exception as e:
                    print(e)

    def test_TCR_geom_methods(self):
        parser = TCRParser.TCRParser()
        test_file = "./STCRpy/test/test_files/8gvb.cif"
        tcr = list(parser.get_tcr_structure("8gvb", test_file).get_TCRs())[0]
        geometry = tcr.calculate_docking_geometry()
        assert "scanning_angle" in geometry
