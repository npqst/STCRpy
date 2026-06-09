import unittest
import glob

import stcrpy
import stcrpy.tcr_processing
from stcrpy.tcr_formats import tcr_haddock

PARSER_FILES = "./test_files/TCRParser_test_files"
HADDOCK_OUT = "./test_files/out/haddock/"


class TestTCRFormats(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tcrs = (
            stcrpy.load_TCRs(glob.glob(f"{PARSER_FILES}/*.pdb"))
            + stcrpy.load_TCRs(glob.glob(f"{PARSER_FILES}/*.cif"))
        )

    def test_tcr_to_haddock(self):
        haddock_formatter = tcr_haddock.HADDOCKFormatter(HADDOCK_OUT)
        for tcr in self.tcrs:
            if tcr is not None:
                haddock_formatter.tcr_to_haddock(tcr)

    def test_pMHC_to_haddock(self):
        haddock_formatter = tcr_haddock.HADDOCKFormatter(save_dir=HADDOCK_OUT)
        for tcr in self.tcrs:
            if (
                tcr is not None
                and len(tcr.get_MHC()) > 0
                and len(tcr.antigen) > 0
                and isinstance(tcr.antigen[0], stcrpy.tcr_processing.AGchain.AGchain)
            ):
                haddock_formatter.pMHC_to_haddock(tcr.get_MHC()[0], tcr.antigen)

    # The two tests below require the full HADDOCK output directory
    # (test_files/TCRHaddock_test_files/387937-tcr_6eqa_mel5_bulged/) which is
    # excluded from git due to size. Re-enable once that data is available.

    # def test_from_haddock_to_TCR_pMHC(self):
    #     from stcrpy.tcr_formats import tcr_haddock
    #     haddock_results_parser = tcr_haddock.HADDOCKResultsParser(
    #         haddock_results_dir="./test_files/TCRHaddock_test_files/387937-tcr_6eqa_mel5_bulged",
    #         tcr_renumbering_file="./test_files/TCRHaddock_test_files/6eqa_TCR_haddock_renumbering.txt",
    #         pmhc_renumbering_file="./test_files/TCRHaddock_test_files/6eqa_pMHC_haddock_renumbering.txt",
    #     )
    #     haddock_results_parser.renumber_all_haddock_predictions()
    #     for file_path in glob.glob(
    #         "./test_files/TCRHaddock_test_files/387937-tcr_6eqa_mel5_bulged/structures/it1/complex*.pdb"
    #     ):
    #         renumbered = file_path.replace("complex", "renumbered_complex")
    #         assert os.path.exists(renumbered)

    # def test_get_haddock_scores(self):
    #     from stcrpy.tcr_formats import tcr_haddock
    #     haddock_results_parser = tcr_haddock.HADDOCKResultsParser(
    #         haddock_results_dir="./test_files/TCRHaddock_test_files/387937-tcr_6eqa_mel5_bulged",
    #     )
    #     scores = haddock_results_parser.get_haddock_scores()
    #     assert len(scores) == 200
