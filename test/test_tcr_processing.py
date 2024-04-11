import unittest
from ..TCRpy.tcr_processing import TCRParser


class TestTCRParser(unittest.TestCase):
    def test_get_tcr_structure_class_I(self):
        parser = TCRParser.TCRParser()

        pdb_file = 'TCRpy/test/test_files/5hyj.pdb'
        tcr = parser.get_tcr_structure('test', pdb_file)
        assert set([''.join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(['DE', 'IJ'])
        assert set([''.join(sorted(x.id)) for x in tcr.get_MHCs()]) == set(['FG', 'AB'])
        assert set([''.join(sorted(x.id)) for x in tcr.get_antigens()]) == set(['C', 'H'])

    def test_get_tcr_structure_class_II(self):
        parser = TCRParser.TCRParser()

        pdb_file = 'TCRpy/test/test_files/6r0e.cif'
        tcr = parser.get_tcr_structure('test', pdb_file)
        assert set([''.join(sorted(x.id)) for x in tcr.get_TCRs()]) == set(['DE'])
        assert set([''.join(sorted(x.id)) for x in tcr.get_MHCs()]) == set(['AB'])
        assert set([''.join(sorted(x.id)) for x in tcr.get_antigens()]) == set(['C'])
