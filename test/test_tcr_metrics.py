import unittest
import os
import glob

import pandas as pd

import stcrpy
import stcrpy.tcr_metrics
from stcrpy.tcr_geometry.TCRAngle import TCRAngle

STRUCTURES = "./test_files/structures"
RMSD_FILES = "./test_files/TCRRMSD_test_files"
HADDOCK_FILES = "./test_files/TCRHaddock_test_files"
IRMSD_FILES = "./test_files/TCRInterfaceRMSD_test_files"


class TestTCRMetrics(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ab_tcr = stcrpy.load_TCRs(f"{STRUCTURES}/8gvb.cif")[0]
        cls.gd_tcr = stcrpy.load_TCRs(f"{STRUCTURES}/8jbv.cif")[0]

    def test_tcr_rmsd(self):
        true_tcr, pred_tcr = stcrpy.load_TCRs(
            {
                "true_7su9": f"{RMSD_FILES}/true_7su9_0_ED.pdb",
                "pred_7su9": f"{RMSD_FILES}/pred_7su9.pdb",
            }
        ).values()

        from stcrpy.tcr_metrics import RMSD

        rmsds = RMSD().calculate_rmsd(pred_tcr, true_tcr, save_alignment=True)

        correct_rmsd = {
            "B": 0.4232598777149801,
            "FWB": 0.4059422,
            "CDRB1": 0.36291695,
            "CDRB2": 0.40915382,
            "CDRB3": 0.5511935,
            "A": 0.6546116886811664,
            "FWA": 0.692916,
            "CDRA1": 0.45867094,
            "CDRA2": 0.4696773,
            "CDRA3": 0.5279137,
        }
        assert all([abs(correct_rmsd[k] - rmsds[k]) < 0.00001 for k in rmsds])

    def test_tcr_rmsd_from_file_list(self):
        target_file_path = f"{RMSD_FILES}/true_structures"
        prediction_file_path = f"{RMSD_FILES}/pred_structures/"

        target_files = sorted(os.listdir(target_file_path))
        prediction_files = sorted(os.listdir(prediction_file_path))

        files = list(
            zip(
                [os.path.join(prediction_file_path, f) for f in prediction_files],
                [os.path.join(target_file_path, f) for f in target_files],
            )
        )

        from stcrpy.tcr_metrics import RMSD

        rmsd_df = RMSD().rmsd_from_files(files)
        assert len(rmsd_df) == 3

        correct_rmsd = pd.read_csv(
            f"{RMSD_FILES}/rmsd_testing.csv",
            index_col="Unnamed: 0",
        )
        map_column_names = lambda c: c.lower() if len(c) > 1 else c
        for idx, rmsd_row in rmsd_df.iterrows():
            reference_row = correct_rmsd[correct_rmsd.pdb == idx]
            if len(reference_row) == 0:
                assert idx == "7sg0"
                continue
            assert all(
                abs(reference_row[map_column_names(col)].item() - rmsd_row[col]) < 0.00001
                for col in rmsd_row.index
            )

    def test_interface_rmsd(self):
        from stcrpy.tcr_metrics import InterfaceRMSD

        dock_files = sorted(
            glob.glob(f"{HADDOCK_FILES}/predictions/renumbered_complex_*.pdb")
        )
        docked_tcrs = stcrpy.load_TCRs(dock_files)
        reference_tcr = stcrpy.load_TCRs(f"{IRMSD_FILES}/6eqa.cif")[0]

        irmsds = [
            InterfaceRMSD().get_interface_rmsd(tcr, reference_tcr) for tcr in docked_tcrs
        ]
        assert len(irmsds) == 4
        detached_peptide_indices = [13, 35, 37, 127]
        assert all(
            r > 0.0
            for i, r in enumerate(irmsds)
            if not any(
                f"renumbered_complex_{p_idx}" in dock_files[i]
                for p_idx in detached_peptide_indices
            )
        )

    def test_dockq(self):
        from stcrpy.tcr_metrics.tcr_dockq import TCRDockQ

        dock_files = sorted(
            glob.glob(f"{HADDOCK_FILES}/predictions/renumbered_complex_*.pdb")
        )[:2]
        docked_tcrs = stcrpy.load_TCRs(dock_files)
        reference_tcr = stcrpy.load_TCRs(f"{IRMSD_FILES}/6eqa.cif")[0]

        dockq_results = [TCRDockQ().tcr_dockq(tcr, reference_tcr) for tcr in docked_tcrs]

        self.assertAlmostEqual(dockq_results[0]['best_result']['TM']['DockQ'], 0.825, places=3)
        self.assertAlmostEqual(dockq_results[0]['best_result']['TM']['iRMSD'], 0.827, places=3)
        self.assertAlmostEqual(dockq_results[0]['best_result']['TM']['LRMSD'], 4.025, places=3)
        self.assertAlmostEqual(dockq_results[0]['best_result']['TM']['fnat'], 0.892, places=3)
        self.assertAlmostEqual(dockq_results[0]['best_result']['TM']['fnonnat'], 0.266, places=3)
        self.assertAlmostEqual(dockq_results[0]['best_result']['TM']['F1'], 0.806, places=3)
        self.assertAlmostEqual(dockq_results[0]['best_result']['TM']['clashes'], 0, places=3)

    def test_TCRAngles(self):
        tcr_angle = TCRAngle()

        angles = tcr_angle.calculate_angles(self.ab_tcr)
        assert self.ab_tcr.get_TCR_angles() == angles

        angles = tcr_angle.calculate_angles(self.gd_tcr)
        assert self.gd_tcr.get_TCR_angles() == angles
