import unittest

import stcrpy

try:
    from stcrpy.tcr_datasets.tcr_graph_dataset import TCRGraphConstructor
    HAS_ML_PKGS = True
except ImportError:
    HAS_ML_PKGS = False

STRUCTURES = "./test_files/structures"


@unittest.skipUnless(HAS_ML_PKGS, "ml_datasets extras not installed (torch/torch_geometric)")
class TestTCRDatasets(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tcr = stcrpy.load_TCRs(f"{STRUCTURES}/8gvb.cif")[0]

    def test_TCRGraphConstructor(self):
        graph_constructor = TCRGraphConstructor()
        assert graph_constructor.config == {
            "node_level": "residue",
            "residue_coord": ["CA"],
            "node_features": "one_hot",
            "edge_features": "distance",
            "tcr_regions": ["all"],
            "include_antigen": True,
            "include_mhc": True,
            "mhc_distance_threshold": 15.0,
        }
        graph_constructor.build_graph(self.tcr)