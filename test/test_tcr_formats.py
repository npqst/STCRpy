import unittest
import glob
import os

import stcrpy.tcr_processing


class TestTCRFormats(unittest.TestCase):
    def test_tcr_to_haddock(self):
        import stcrpy

        tcrs = stcrpy.load_TCRs(glob.glob("./test_files/*.pdb")) + stcrpy.load_TCRs(
            glob.glob("./test_files/*.cif")
        )

        from stcrpy.tcr_formats import tcr_haddock

        haddock_formatter = tcr_haddock.HADDOCKFormatter("./test_files/out/haddock/")
        for tcr in tcrs:
            haddock_formatter.tcr_to_haddock(tcr)

    def test_pMHC_to_haddock(self):
        import stcrpy

        tcrs = stcrpy.load_TCRs(glob.glob("./test_files/*.pdb")) + stcrpy.load_TCRs(
            glob.glob("./test_files/*.cif")
        )

        from stcrpy.tcr_formats import tcr_haddock

        haddock_formatter = tcr_haddock.HADDOCKFormatter(
            save_dir="./test_files/out/haddock/"
        )
        for tcr in tcrs:
            if (
                len(tcr.get_MHC()) > 0
                and len(tcr.antigen) > 0
                and isinstance(tcr.antigen[0], stcrpy.tcr_processing.AGchain.AGchain)
            ):
                haddock_formatter.pMHC_to_haddock(tcr.get_MHC()[0], tcr.antigen)

    def test_from_haddock_to_TCR_pMHC(self):
        import stcrpy
        from stcrpy.tcr_formats import tcr_haddock

        haddock_results_parser = tcr_haddock.HADDOCKResultsParser(
            haddock_results_dir="./test_files/TCRHaddock_test_files/387937-tcr_6eqa_mel5_bulged",
            tcr_renumbering_file="./test_files/TCRHaddock_test_files/6eqa_TCR_haddock_renumbering.txt",
            pmhc_renumbering_file="./test_files/TCRHaddock_test_files/6eqa_pMHC_haddock_renumbering.txt",
        )

        for renumbered_file_path in glob.glob(
            "./test_files/TCRHaddock_test_files/387937-tcr_6eqa_mel5_bulged/structures/it1/renumbered_complex*.pdb"
        ):
            os.remove(renumbered_file_path)

        haddock_results_parser.renumber_all_haddock_predictions()
        for file_path in glob.glob(
            "./test_files/TCRHaddock_test_files/387937-tcr_6eqa_mel5_bulged/structures/it1/complex*.pdb"
        ):
            renumbered_file_path = (
                file_path.split("complex")[0]
                + "renumbered_complex"
                + file_path.split("complex")[1]
            )
            assert os.path.exists(renumbered_file_path)

    def test_get_haddock_scores(self):
        import stcrpy
        from stcrpy.tcr_formats import tcr_haddock

        haddock_results_parser = tcr_haddock.HADDOCKResultsParser(
            haddock_results_dir="./test_files/TCRHaddock_test_files/387937-tcr_6eqa_mel5_bulged",
        )

        scores = haddock_results_parser.get_haddock_scores()
        assert len(scores) == 200
