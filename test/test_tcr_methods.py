import unittest

import stcrpy
from stcrpy import fetch_TCRs

STRUCTURES = "./test_files/structures"


class TestTCRMethods(unittest.TestCase):
    def test_fetch_tcr(self):
        # Tests the fetch_TCRs function end-to-end: download, parse, and warn on non-TCR
        tcrs = fetch_TCRs("6eqa")
        self.assertIsInstance(tcrs[0], stcrpy.tcr_processing.abTCR)

        with self.assertWarns(UserWarning):
            non_tcr = fetch_TCRs("8zt4")
        self.assertEqual(non_tcr, [])

    def test_load_tcr(self):
        tcrs = stcrpy.load_TCRs(f"{STRUCTURES}/6eqa.cif")
        self.assertIsInstance(tcrs[0], stcrpy.tcr_processing.abTCR)