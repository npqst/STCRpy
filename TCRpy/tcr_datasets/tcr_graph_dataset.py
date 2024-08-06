
import warnings
import numpy as np
import itertools
import os
import os.path as osp
import sys
from typing import List, Any
import re

from torch_geometric.data import Data, DataLoader, Dataset
from torch_geometric.io import fs
import torch
import torch.nn as nn 
import torch.nn.functional as F

from ..tcr_processing import TCR, TCRParser
from . import utils


class TCRGraphConstructor():
    def __init__(self, config=None):
        if config is None:
            config = {
                'node_level': 'residue',
                'residue_coord': ['CA'],
                'node_features': 'one_hot',
                'edge_features': 'distance',
                'tcr_regions': ['all'],
                'include_antigen': True,
                'include_mhc': True,
                'mhc_distance_threshold': 15.,
            }

        # assert that minimum amount of configuration is set
        assert len(set(['node_level', 'node_features', 'edge_features']) - set(config.keys())) == 0

        self.config = config

        self.node_selector = self._get_node_selector()
        self.node_featuriser = self._get_node_featuriser()
        self.edge_featuriser = self._get_edge_featuriser()

    def set_node_selector(self, node_selector_function):
        import Bio
        test_res = Bio.PDB.Residue.Residue(id=(' ', 3, ' '), resname='GLY', segid=' ')
        atom_N = Bio.PDB.Atom.Atom(name='N', coord=np.array([23.399, -5.842, 19.395]), bfactor=67.02)
        atom_O = Bio.PDB.Atom.Atom(name='O', coord=np.array([24.17, -8.195, 21.998]), bfactor=67.02)
        atom_C = Bio.PDB.Atom.Atom(name='C', coord=np.array([23.617, -7.263, 21.414]), bfactor=67.02)
        atom_CA = Bio.PDB.Atom.Atom(name='CA', coord=np.array([24.316, -6.528, 20.288]), bfactor=67.02)
        for a in [atom_N, atom_O, atom_C, atom_CA]:
            test_res.add(a)
        try:
            node_selector_function(test_res)
        except Exception as e:
            raise ValueError(f'Node selector function should generate node from Bio.PDB.Residue instance. Raised error {e}')
        self.node_selector = node_selector_function

    def set_node_featuriser(self, node_featuriser_function, test_input=None):
        if test_input is None:
            warnings.warn('No test input provided for new node featuriser, using Bio.PDB.Atom instance')
            import Bio
            test_input = Bio.PDB.Atom.Atom(name='CA', coord=np.array([24.316, -6.528, 20.288]), bfactor=67.02)
        try:
            features = node_featuriser_function(test_input)
        except Exception as e:
            raise ValueError(f'Node featuriser function could not featurise node {test_input}. Raised error {e}')
        assert isinstance(features, torch.tensor) or features is None, 'Node featuriser should generate torch tensor'
        self.node_featuriser = node_featuriser_function

    def set_edge_featuriser(self, edge_featuriser_function, test_input=None):
        if test_input is None:
            warnings.warn('No test input provided for new edge featuriser, using Bio.PDB.Atom instance')
            import Bio
            test_input = [
                Bio.PDB.Atom.Atom(name='CA', coord=np.array([24.316, -6.528, 20.288]), bfactor=67.02),
                Bio.PDB.Atom.Atom(name='CA', coord=np.array([27.623, -12.28, 23.288]), bfactor=67.02),
                Bio.PDB.Atom.Atom(name='CA', coord=np.array([16.36, 8.58, 30.288]), bfactor=67.02)
                ]
        try:
            edges, edge_features, edge_weights = edge_featuriser_function(test_input)
        except Exception as e:
            raise ValueError(f'Edge featuriser function could not featurise edge {test_input}. Raised error {e}')
        assert edges.shape[0] == 2, f'Edge indices must be 2D to define connected nodes. Edge shape was {edges.shape}'
        if edge_features is not None:
            assert edges.shape[1] == len(edge_features) and isinstance(edge_features, torch.tensor), 'Edge features configuration invalid'
        if edge_weights is not None:
            assert edges.shape[1] == len(edge_weights) and isinstance(edge_weights, torch.tensor), 'Edge weights configuration invalid'
        
        self.edge_featuriser = edge_featuriser_function


    def _calculate_distance_matrix(self, coord_1, coord_2):
        assert coord_1.shape[-1] == coord_2.shape[-1] == 3
        coord_1_matrix = np.tile(coord_1, (len(coord_2), 1, 1))
        coord_2_matrix = np.moveaxis(np.tile(coord_2, (len(coord_1), 1, 1)), 0, 1)
        assert coord_1_matrix.shape == coord_2_matrix.shape

        euclidian_dist_mat = np.sqrt(np.sum((coord_1_matrix - coord_2_matrix) ** 2, axis=-1))
        return euclidian_dist_mat.squeeze()

    def _get_node_selector(self):
        if self.config['node_level'] == 'residue':
            if 'residue_coord' not in self.config or self.config['residue_coord'] == ['CA']:
                # generate single node per residue with coordinate of CA atom.
                def node_selector(residue):
                    if 'CA' in residue.child_dict:
                        return [residue['CA']]
                    else:
                        return [None]
                return node_selector
        else:
            NotImplementedError

    def _get_node_featuriser(self):
        if self.config['node_features'] == 'one_hot':
            def one_hot_encoding(atom_node):
                """one_hot_encoding consists of ...
                4 dims for chain type: [TCR alpha, TCR beta, peptide, MHC]
                7 dims for CDR loop: [Not CDR, CDRA1, CDRA2, CDRA3, CDRB1, CDRB2, CDRB2]
                20 dims for residue encoding: ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'XXX'
                37 dims for atom encoding
                Args:
                    atom_node (_type_): _description_
                """
                if hasattr(atom_node.parent.parent, 'chain_type'):
                    chain_type = atom_node.parent.parent.chain_type
                else:
                    if hasattr(atom_node.parent.parent, 'type'): 
                        chain_type = atom_node.parent.parent.type
                        if chain_type == 'peptide':
                            chain_type = 'antigen'
                    else:
                        chain_type = atom_node.parent.parent.MHC_type
                        chain_type = 'MHC'      # if calling MHC_type doesn't raise an error, chain is MHC
                chain_type = 'MHC' if chain_type in utils.MHC_CHAIN_TYPES else chain_type

                if chain_type in utils.CHAIN_TYPE_ONE_HOT_ENCODING:
                    chain_type_onehot_encoding = utils.CHAIN_TYPE_ONE_HOT_ENCODING[chain_type]     # one hot as integer
                else: 
                    warnings.warn(f'Could not resolve chain type: {chain_type} of node: {atom_node}|{atom_node.parent.parent}|{atom_node.parent.parent.parent.id}')
                    return None

                if chain_type in ['A', 'B', 'G', 'D']:
                    region = atom_node.parent.region.capitalize()
                else:
                    region = 'NOT_CDR'
                if region in utils.TCR_REGION_ONE_HOT_ENCODING:
                    region_onehot_encoding = utils.TCR_REGION_ONE_HOT_ENCODING[region]
                else:
                    region_onehot_encoding = utils.TCR_REGION_ONE_HOT_ENCODING['NOT_CDR']

                residue_onehot_encoding = utils.AMINO_ACID_ONEHOT_ENCODING[atom_node.parent.resname.strip()]

                atom37_onehot_encoding = utils.ATOM37_ATOM_ONEHOT_ENCODING[atom_node.fullname.strip()]
                atom_onehot_encoding = torch.concat([
                    F.one_hot(torch.tensor(chain_type_onehot_encoding), num_classes=4),
                    F.one_hot(torch.tensor(region_onehot_encoding), num_classes=7),
                    F.one_hot(torch.tensor(residue_onehot_encoding), num_classes=21),
                    F.one_hot(torch.tensor(atom37_onehot_encoding), num_classes=37),
                ])
                return atom_onehot_encoding
            return one_hot_encoding
        else:
            raise NotImplementedError('Node featurisation method not recognised')

    def _get_edge_featuriser(self):
        if self.config['edge_features'] == 'distance':
            import scipy

            def distance_edges(nodes, distance_cutoff=15.):
                dist_mat = np.triu(np.zeros((len(nodes), len(nodes))))
                coords = np.asarray([a.get_coord() for a in nodes])
                dist_mat[np.arange(len(nodes))[:, None] < np.arange(len(nodes))] = scipy.spatial.distance.pdist(coords)
                dist_mat = dist_mat + dist_mat.T + (distance_cutoff * np.eye(len(nodes)))      # add diagonal to remove self edges
                edges = np.argwhere(dist_mat < distance_cutoff)
                edge_features = dist_mat[edges[:, 0], edges[:, 1]]
                return torch.from_numpy(edges), torch.from_numpy(edge_features), None
            
            return distance_edges

        else: 
            raise NotImplementedError('Edge featurisation method not recognised')

    def build_graph(self, tcr: TCR, label=None):
        nodes = []
        coordinates = []
        
        
        if 'tcr_regions' not in self.config or self.config['tcr_regions'] == ['all'] or self.config['tcr_regions'] is None:
            tcr_nodes = [a for res in tcr.get_residues() for a in self.node_selector(res) if res.id[0].strip() == '']     # filters out waters and other
            tcr_coords = np.array([a.get_coord() for a in tcr_nodes])
            nodes.extend(tcr_nodes)
            coordinates.extend(tcr_coords)

        if 'include_antigen' in self.config and self.config['include_antigen']:
            if len(tcr.get_antigen()) == 0:
                warnings.warn(f'No antigen found for TCR {tcr}. Antigen not included in graph.')
            else:
                antigen_nodes = [a for res in tcr.get_antigen()[0].get_residues() for a in self.node_selector(res) if res.id[0].strip() == '']
                antigen_coords = np.array([a.get_coord() for a in antigen_nodes])
                nodes.extend(antigen_nodes)
                coordinates.extend(antigen_coords)
        
        if 'include_mhc' in self.config and self.config['include_mhc']:
            if len(tcr.get_MHC()) == 0:
                warnings.warn(f'No MHC found for TCR {tcr}. MHC not included in graph.')
            else:
                mhc_nodes = [a for res in tcr.get_MHC()[0].get_residues() for a in self.node_selector(res) if res.id[0].strip() == '']
                mhc_coords = np.array([a.get_coord() for a in mhc_nodes])
                if 'mhc_distance_threshold' in self.config:
                    dist_mat = self._calculate_distance_matrix(tcr_coords, mhc_coords)
                    mhc_node_mask = np.sum(dist_mat < self.config['mhc_distance_threshold'], axis=-1) > 0           # shape is (len(mhc_nodes), len(tcr_nodes))
                    mhc_nodes = list(itertools.compress(mhc_nodes, mhc_node_mask))
                    mhc_coords = np.array(list(itertools.compress(mhc_coords, mhc_node_mask)))
                nodes.extend(mhc_nodes)
                coordinates.extend(mhc_coords)
        
        node_features = [self.node_featuriser(n) for n in nodes]
        
        # remove nodes that could not be featurised
        indices_to_remove = [idx for idx, n in enumerate(node_features) if n is None]        
        for idx in indices_to_remove:
            nodes.pop(idx)
            node_features.pop(idx)
        if len(indices_to_remove) > 0:
            warnings.warn(f'{len(indices_to_remove)} nodes removed from original node list of TCR {tcr.parent.parent.id}_{tcr.id}')
        assert len(nodes) == len(node_features)
        
        edge_index, edge_features, edge_weight = self.edge_featuriser(nodes)

        graph = Data(
            x=node_features,
            edge_index=edge_index.T,
            edge_attr=edge_features,
            edge_weight=edge_weight,
            pos=torch.from_numpy(np.array(coordinates)),
            y=label,
            tcr_id=f'{tcr.parent.parent.id}_{tcr.id}',
        )
        return graph


class TCRGraphDataset(Dataset):
    def __init__(self, root=None, graph_config=None, *args, **kwargs):
        
        self.graph_constructor = TCRGraphConstructor(config=graph_config)

        if isinstance(root, str):
            self.root = root
            self._raw_file_names = [
                os.path.join(root, 'raw_files', f)
                for f in os.listdir(os.path.join(root, 'raw_files'))
                if f.endswith('.pdb') or f.endswith('.cif')
                ]
            
        self._processed_file_names = []
        
        super(TCRGraphDataset, self).__init__(root=root, *args, **kwargs)

    @property
    def raw_file_names(self):
        return self._raw_file_names

    @property
    def processed_file_names(self):
        return self._processed_file_names

    @staticmethod
    def _tcr_generator(tcr_parser, tcr_pdb_iter):
        for tcr in tcr_pdb_iter:
            tcr_id = tcr.split('/')[-1].split('.')[0]
            yield tcr_parser.get_tcr_structure(tcr_id, tcr).get_TCRs()

    def process(self):
        tcr_parser = TCRParser.TCRParser()
        for tcr_object in self._tcr_generator(tcr_parser, self.raw_file_names):
            for tcr in tcr_object:
                try:
                    tcr_graph = self.graph_constructor.build_graph(tcr)
                    processed_file_path = os.path.join(self.root, 'processed', f'{tcr_graph.tcr_id}.pt')
                    
                    torch.save(tcr_graph, processed_file_path)
                    self._processed_file_names.append(processed_file_path)
                except Exception as e:
                    warnings.warn(f'Dataset parsing error: {str(e)} for TCR: {tcr}')
    
    def len(self):
        return len(self._processed_file_names)
    
    def get(self, idx):
        graph = torch.load(self._processed_file_names[idx])
        return graph
