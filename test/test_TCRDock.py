import unittest
from ..src.tcr_processing import TCRParser
from ..src.tcr_geometry import TCRDock


class TestTCRDock(unittest.TestCase):
    def test_TCRDock_init(self):
        parser = TCRParser.TCRParser()
        pdb_file = 'TCRpy/test/test_files/5hyj.pdb'
        tcr = parser.get_tcr_structure('test', pdb_file)
        [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]

    def test_calculate_docking_angle(self):
        parser = TCRParser.TCRParser()
        pdb_file = 'TCRpy/test/test_files/5hyj.pdb'
        tcr = parser.get_tcr_structure('test', pdb_file)

        tcr_docks = [TCRDock.TCRDock(x) for x in tcr.get_TCRs()]
        assert all(
            [
                x.calculate_docking_angle() < 50.
                and x.calculate_docking_angle() > 40.
                for x in tcr_docks
            ])
