import unittest

import os
import pandas as pd
import glob


class TestTCRMethods(unittest.TestCase):
    def test_fetch_tcr(self):
        import stcrpy
        from stcrpy import fetch_TCR

        tcr = fetch_TCR("6eqa")
        self.assertIsInstance(tcr, stcrpy.tcr_processing.abTCR)

        with self.assertWarns(UserWarning):
            non_tcr = fetch_TCR("8zt4")
        self.assertIsNone(non_tcr)
