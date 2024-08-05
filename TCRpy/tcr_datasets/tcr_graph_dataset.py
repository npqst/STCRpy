
import warnings
import numpy as np
import itertools

from torch_geometric.data import Data, DataLoader
import torch
import torch.nn as nn 
import torch.nn.functional as F

from ..tcr_processing import TCR
from . import utils


class TCRGraph():
    def __init__(self, tcr: TCR, config:dict=None):
        self._nodes = []
        self._node_coords = []
        self._edge_ind = []
        self._edge_attr = []
        self._id = tcr.id()

        self._graph_constructor = TCRGraphConstructor(config=config)

        self._construct_graph(tcr, config=config)
    
    def _construct_graph(self, tcr:TCR):
        (
            self._nodes,
            self._node_coords, 
            self._edge_ind,
            self._edge_attr,
        ) = self._graph_constructor.build_graph(tcr)
        


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
                    return [residue['CA']]
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
                chain_type_onehot_encoding = utils.CHAIN_TYPE_ONE_HOT_ENCODING[chain_type]     # one hot as integer
                if chain_type in ['A', 'B', 'G', 'D']:
                    region = atom_node.parent.region.capitalize()
                else:
                    region = 'NOT_CDR'
                if region in utils.TCR_REGION_ONE_HOT_ENCODING:
                    region_onehot_encoding = utils.TCR_REGION_ONE_HOT_ENCODING[region]
                else:
                    region_onehot_encoding = utils.TCR_REGION_ONE_HOT_ENCODING['NOT_CDR']

                residue_onehot_encoding = utils.AMINO_ACID_ONEHOT_ENCODING[atom_node.parent.resname]

                atom37_onehot_encoding = utils.ATOM37_ATOM_ONEHOT_ENCODING[atom_node.fullname]
                atom_onehot_encoding = torch.concat([
                    F.one_hot(torch.tensor(chain_type_onehot_encoding), num_classes=4),
                    F.one_hot(torch.tensor(region_onehot_encoding), num_classes=7),
                    F.one_hot(torch.tensor(residue_onehot_encoding), num_classes=21),
                    F.one_hot(torch.tensor(atom37_onehot_encoding), num_classes=37),
                ])
                return atom_onehot_encoding
            return one_hot_encoding

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
                return torch.from_numpy(edges), torch.from_numpy(edge_features)
            
            return distance_edges

    def build_graph(self, tcr: TCR, label=None):
        nodes = []
        coordinates = []
        node_selector = self._get_node_selector()
        
        if 'tcr_regions' not in self.config or self.config['tcr_regions'] == ['all'] or self.config['tcr_regions'] is None:
            tcr_nodes = [a for res in tcr.get_residues() for a in node_selector(res)]
            tcr_coords = np.array([a.get_coord() for a in tcr_nodes])
            nodes.extend(tcr_nodes)
            coordinates.extend(tcr_coords)

        if 'include_antigen' in self.config and self.config['include_antigen']:
            if len(tcr.get_antigen()) == 0:
                warnings.warn(f'No antigen found for TCR {tcr}. Antigen not included in graph.')
            else:
                antigen_nodes = [a for res in tcr.get_antigen()[0].get_residues() for a in node_selector(res) ]
                antigen_coords = np.array([a.get_coord() for a in antigen_nodes])
                nodes.extend(antigen_nodes)
                coordinates.extend(antigen_coords)
        
        if 'include_mhc' in self.config and self.config['include_mhc']:
            if len(tcr.get_MHC()) == 0:
                warnings.warn(f'No MHC found for TCR {tcr}. MHC not included in graph.')
            else:
                mhc_nodes = [a for res in tcr.get_MHC()[0].get_residues() for a in node_selector(res)]
                mhc_coords = np.array([a.get_coord() for a in mhc_nodes])
                if 'mhc_distance_threshold' in self.config:
                    dist_mat = self._calculate_distance_matrix(tcr_coords, mhc_coords)
                    mhc_node_mask = np.sum(dist_mat < self.config['mhc_distance_threshold'], axis=-1) > 0           # shape is (len(mhc_nodes), len(tcr_nodes))
                    mhc_nodes = list(itertools.compress(mhc_nodes, mhc_node_mask))
                    mhc_coords = np.array(list(itertools.compress(mhc_coords, mhc_node_mask)))
                nodes.extend(mhc_nodes)
                coordinates.extend(mhc_coords)
        
        node_featuriser = self._get_node_featuriser()
        node_features = [node_featuriser(n) for n in nodes]
        
        edge_featuriser = self._get_edge_featuriser()
        edge_index, edge_features = edge_featuriser(nodes)

        graph = Data(
            x=node_features,
            edge_index=edge_index.T,
            edge_attr=edge_features,
            pos=torch.from_numpy(np.array(coordinates)),
            y=label,
            tcr_id=tcr.id,
        )
        return graph

