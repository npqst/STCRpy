import unittest

import os
import stcrpy.tcr_methods
import pandas as pd
import glob


import stcrpy
import stcrpy.tcr_metrics


class TestTCRMetrics(unittest.TestCase):
    def test_tcr_rmsd(self):
        true_tcr, pred_tcr = stcrpy.load_TCRs(
            {
                "true_7su9": "./test_files/TCRRMSD_test_files/true_7su9_0_ED.pdb",
                "pred_7su9": "./test_files/TCRRMSD_test_files/pred_7su9.pdb",
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
        target_file_path = "./test_files/TCRRMSD_test_files/true_structures"
        target_files = os.listdir(target_file_path)
        prediction_file_path = "./test_files/TCRRMSD_test_files/pred_structures/"
        prediction_files = os.listdir(prediction_file_path)

        target_files.sort()
        prediction_files.sort()

        files = list(
            zip(
                [os.path.join(prediction_file_path, f) for f in prediction_files],
                [os.path.join(target_file_path, f) for f in target_files],
            )
        )

        from stcrpy.tcr_metrics import RMSD

        rmsd_df = RMSD().rmsd_from_files(files)
        assert len(rmsd_df) == 46

        correct_rmsd = pd.read_csv(
            "test_files/TCRRMSD_test_files/rmsd_testing.csv",
            index_col="Unnamed: 0",
        )
        map_column_names = lambda c: c.lower() if len(c) > 1 else c
        for idx, rmsd_row in rmsd_df.iterrows():
            reference_row = correct_rmsd[correct_rmsd.pdb == idx]
            if len(reference_row) == 0:
                assert idx == "7sg0"
                continue
            assert all(
                [
                    abs(reference_row[map_column_names(col)].item() - rmsd_row[col])
                    < 0.00001
                    for col in rmsd_row.index
                ]
            )

    def test_interface_rmsd(self):
        from stcrpy.tcr_metrics import InterfaceRMSD

        interface_rmsd = InterfaceRMSD()

        dock_files = sorted(
            glob.glob(
                "./test_files/TCRHaddock_test_files/387937-tcr_6eqa_mel5_bulged/structures/it1/renumbered_complex_*.pdb"
            )
        )
        docked_tcrs = stcrpy.load_TCRs(dock_files)

        reference_tcr = stcrpy.load_TCRs(
            "./test_files/TCRInterfaceRMSD_test_files/6eqa.cif"
        )[0]

        irmsds = [
            interface_rmsd.get_interface_rmsd(tcr, reference_tcr) for tcr in docked_tcrs
        ]
        assert len(irmsds) == 200
        detached_peptide_indices = [13, 35, 37, 127]
        assert all(
            [
                r > 0.0
                for i, r in enumerate(irmsds)
                if not any(
                    [
                        f"renumbered_complex_{p_idx}" in dock_files[i]
                        for p_idx in detached_peptide_indices
                    ]
                )
            ]
        )
