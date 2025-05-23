import unittest

import os
import pandas as pd
import glob


class TestTCRMethods(unittest.TestCase):
    def test_fetch_tcr(self):
        import stcrpy
        from stcrpy import fetch_TCRs

        tcrs = fetch_TCRs("6eqa")
        self.assertIsInstance(tcr[0], stcrpy.tcr_processing.abTCR)

        with self.assertWarns(UserWarning):
            non_tcr = fetch_TCRs("8zt4")
        self.assertEqual(non_tcr, [])
