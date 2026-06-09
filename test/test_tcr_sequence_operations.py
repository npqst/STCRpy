import unittest

import stcrpy

STRUCTURES = "./test_files/structures"


class TestTCRSequenceOperations(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tcr = stcrpy.load_TCRs(f"{STRUCTURES}/8gvb.cif")[0]

    def test_get_germlines(self):
        germline_info = self.tcr.get_germline_assignments()

        self.assertIsInstance(germline_info, dict)
        self.assertGreater(len(germline_info), 0)
        for chain_id, assignment in germline_info.items():
            self.assertIsInstance(chain_id, str)
            self.assertIsNotNone(assignment)

    def test_get_mhc_alleles(self):
        mhc_info = self.tcr.get_MHC_allele_assignments()

        self.assertIsInstance(mhc_info, list)
        self.assertGreater(len(mhc_info), 0)