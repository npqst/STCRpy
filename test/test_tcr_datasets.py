import unittest

from ..TCRpy.tcr_processing.TCRParser import TCRParser
from ..TCRpy.tcr_datasets.tcr_graph_dataset import TCRGraphConstructor


class TestTCRDatasets(unittest.TestCase):    
    def __init__(self, *args, **kwargs):
        super(TestTCRDatasets, self).__init__(*args, **kwargs)

        parser = TCRParser()
        tcr = list(parser.get_tcr_structure('test', 'TCRpy/test/test_files/8gvb.cif').get_TCRs())[0]
        self.tcr = tcr

    def test_TCRGraphConstructor(self):
        graph_constructor = TCRGraphConstructor()
        assert graph_constructor.config == {
                'node_level': 'residue',
                'residue_coord': ['CA'],
                'node_features': 'one_hot',
                'edge_features': 'distance',
                'tcr_regions': ['all'],
                'include_antigen': True,
                'include_mhc': True,
                'mhc_distance_threshold': 15.,
            }
        
        graph_constructor.build_graph(self.tcr)