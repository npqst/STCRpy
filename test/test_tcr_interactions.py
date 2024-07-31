import unittest
import plip
from plip.structure.preparation import PDBComplex

import openbabel
from rdkit import Chem


from ..TCRpy.tcr_processing import TCRParser
from ..TCRpy.tcr_interactions.TCRpMHC_PLIP_Model_Parser import TCRpMHC_PLIP_Model_Parser
from ..TCRpy.tcr_interactions.PLIPParser import PLIPParser

class TestTCRInteractions(unittest.TestCase):

    def test_tcrpmhc_plip_model_parser(self):

        parser = TCRParser.TCRParser()

        model_parser = TCRpMHC_PLIP_Model_Parser()

        test_file = './TCRpy/test/test_files/8gvb.cif'
        tcr = [x for x in parser.get_tcr_structure('tmp', test_file).get_TCRs()][0]

        mol, renumbering, domains = model_parser.parse_tcr_pmhc_complex(tcr)
        assert isinstance(mol, PDBComplex)
        assert set(renumbering.keys()) == set(['A', 'B', 'C', 'D'])
        assert domains == {'VA': 'A', 'VB': 'B'}
        mol.analyze()

        mol = model_parser.parse_tcr_pmhc_complex(tcr, renumber=False)
        assert isinstance(mol, PDBComplex)
        mol.analyze()

    def test_plip_parser(self):
        parser = TCRParser.TCRParser()
        model_parser = TCRpMHC_PLIP_Model_Parser()
        test_file = './TCRpy/test/test_files/8gvb.cif'
        tcr = [x for x in parser.get_tcr_structure('tmp', test_file).get_TCRs()][0]
        mol, renumbering, domains = model_parser.parse_tcr_pmhc_complex(tcr)
        mol.analyze()

        plip_parser = PLIPParser()
        interactions = plip_parser.parse_complex(mol, tcr, renumbering, domains)
        
        assert len(interactions) == 28
        assert len(interactions[interactions.type == 'hbond']) == 12
        assert len(interactions[interactions.type == 'hydrophobic']) == 12
        assert len(interactions[interactions.type == 'pistack']) == 1
        assert len(interactions[interactions.type == 'saltbridge']) == 3

        assert len(interactions[interactions.domain == 'VB']) == 1
        assert interactions[interactions.domain == 'VB'].protein_residue.item() == 'ASP'
        assert interactions[interactions.domain == 'VB'].protein_number.item() == 96

    def test_TCR_plip_methods(self):
        parser = TCRParser.TCRParser()
        test_file = './TCRpy/test/test_files/8gvb.cif'
        tcr = [x for x in parser.get_tcr_structure('tmp', test_file).get_TCRs()][0]

        interactions = tcr.profile_peptide_interactions()
        
        assert len(interactions) == 28
        assert len(interactions[interactions.type == 'hbond']) == 12
        assert len(interactions[interactions.type == 'hydrophobic']) == 12
        assert len(interactions[interactions.type == 'pistack']) == 1
        assert len(interactions[interactions.type == 'saltbridge']) == 3

        assert len(interactions[interactions.domain == 'VB']) == 1
        assert interactions[interactions.domain == 'VB'].protein_residue.item() == 'ASP'
        assert interactions[interactions.domain == 'VB'].protein_number.item() == 96









