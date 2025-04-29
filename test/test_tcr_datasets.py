import unittest

import stcrpy
from stcrpy.tcr_processing.TCRParser import TCRParser
from stcrpy.tcr_datasets.tcr_graph_dataset import TCRGraphConstructor, TCRGraphDataset


class TestTCRDatasets(unittest.TestCase):

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
        tcr = stcrpy.fetch_TCR("8gvb")
        graph_constructor.build_graph(tcr)

    def test_TCRGraphDataset(self):
        dataset = TCRGraphDataset(
            root="./test_files/TCRGraphDataset_test_files",
            data_paths="./test_files/TCRGraphDataset_test_files/raw_files",
            force_reload=True,
        )
        print(dataset)
        for i in range(len(dataset)):
            datapoint = dataset[i]
            print(datapoint)
