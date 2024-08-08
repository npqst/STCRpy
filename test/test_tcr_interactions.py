import unittest
import pathlib
import warnings
import os

from plip.structure.preparation import PDBComplex

from ..TCRpy.tcr_processing import TCRParser
from ..TCRpy.tcr_interactions.TCRpMHC_PLIP_Model_Parser import TCRpMHC_PLIP_Model_Parser
from ..TCRpy.tcr_interactions.PLIPParser import PLIPParser


class TestTCRInteractions(unittest.TestCase):

    def test_tcrpmhc_plip_model_parser(self):

        parser = TCRParser.TCRParser()

        model_parser = TCRpMHC_PLIP_Model_Parser()

        test_file = "./TCRpy/test/test_files/8gvb.cif"
        tcr = [x for x in parser.get_tcr_structure("tmp", test_file).get_TCRs()][0]

        mol, renumbering, domains = model_parser.parse_tcr_pmhc_complex(tcr)
        assert isinstance(mol, PDBComplex)
        assert set(renumbering.keys()) == set(["A", "B", "C", "D"])
        assert domains == {"VA": "A", "VB": "B"}
        mol.analyze()

        mol = model_parser.parse_tcr_pmhc_complex(tcr, renumber=False)
        assert isinstance(mol, PDBComplex)
        mol.analyze()

    def test_plip_parser(self):
        parser = TCRParser.TCRParser()
        model_parser = TCRpMHC_PLIP_Model_Parser()
        test_file = "./TCRpy/test/test_files/8gvb.cif"
        tcr = [x for x in parser.get_tcr_structure("tmp", test_file).get_TCRs()][0]
        mol, renumbering, domains = model_parser.parse_tcr_pmhc_complex(tcr)
        mol.analyze()

        plip_parser = PLIPParser()
        interactions = plip_parser.parse_complex(mol, tcr, renumbering, domains)

        assert len(interactions) == 28
        assert len(interactions[interactions.type == "hbond"]) == 12
        assert len(interactions[interactions.type == "hydrophobic"]) == 12
        assert len(interactions[interactions.type == "pistack"]) == 1
        assert len(interactions[interactions.type == "saltbridge"]) == 3

        assert len(interactions[interactions.domain == "VB"]) == 1
        assert interactions[interactions.domain == "VB"].protein_residue.item() == "ASP"
        assert interactions[interactions.domain == "VB"].protein_number.item() == 96

    def test_TCR_plip_methods(self):
        parser = TCRParser.TCRParser()
        test_file = "./TCRpy/test/test_files/8gvb.cif"
        tcr = [x for x in parser.get_tcr_structure("tmp", test_file).get_TCRs()][0]

        interactions = tcr.profile_peptide_interactions()

        assert len(interactions) == 28
        assert len(interactions[interactions.type == "hbond"]) == 12
        assert len(interactions[interactions.type == "hydrophobic"]) == 12
        assert len(interactions[interactions.type == "pistack"]) == 1
        assert len(interactions[interactions.type == "saltbridge"]) == 3

        assert len(interactions[interactions.domain == "VB"]) == 1
        assert interactions[interactions.domain == "VB"].protein_residue.item() == "ASP"
        assert interactions[interactions.domain == "VB"].protein_number.item() == 96

    def test_pymol_visualisation(self):
        parser = TCRParser.TCRParser()
        model_parser = TCRpMHC_PLIP_Model_Parser()
        test_file = "./TCRpy/test/test_files/8gvb.cif"
        tcr = [x for x in parser.get_tcr_structure("tmp", test_file).get_TCRs()][0]
        mol, renumbering, domains = model_parser.parse_tcr_pmhc_complex(tcr)
        # mol.analyze()

        plip_parser = PLIPParser()

        # test if plip visualisations are generated

        pymol_plip_session_name = "./ATOM1NGLYB123_PROTEIN_UNL_Z_1.pse"
        if os.path.exists(pymol_plip_session_name):
            os.remove(
                pymol_plip_session_name
            )  # remove so test can check if file is generated

        try:  # test that should complete if pymol is installed
            import pymol

            try:
                plip_parser._visualize_interactions(mol)
            except (
                pymol.CmdException
            ):  # sometimes function needs to run twice? Probably due to pymol loading and object selection latency
                plip_parser._visualize_interactions(mol)
            path = pathlib.Path(pymol_plip_session_name)
            assert path.is_file()
            os.remove(pymol_plip_session_name)  # cleans up after test
        except (
            ModuleNotFoundError
        ) as e:  # test that shoudld complete if pymol is not installed
            if "pymol" not in str(e):
                raise ValueError("Only except pymol not found errors")
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                plip_parser._visualize_interactions(mol)
                assert len(w) == 1  # check only one warning raised
                # check warning tells user to install pymol
                assert (
                    "conda install -c conda-forge -c schrodinger pymol-bundle"
                    in str(w[0].message)
                )

    def test_create_pymol_session(self):
        parser = TCRParser.TCRParser()
        test_file = "./TCRpy/test/test_files/8gvb.cif"
        tcr = [x for x in parser.get_tcr_structure("test_8gvb", test_file).get_TCRs()][
            0
        ]

        plip_parser = PLIPParser()

        try:  # tests to run if pymol is installed
            import pymol

            # test automated file name generation
            saved_session = plip_parser.create_pymol_session(tcr)
            assert saved_session == f"{tcr.parent.parent.id}_{tcr.id}_interactions.pse"
            assert pathlib.Path(saved_session).exists()
            os.remove(saved_session)  # clean up after test

            # test saving to specified file
            session_file = "./TCRpy/test/test_files/out/interactions/8gvb_test.pse"
            saved_session = plip_parser.create_pymol_session(tcr, save_as=session_file)
            assert pathlib.Path(saved_session).exists()
            assert session_file == saved_session

            # check tmp file clean up works
            tmp_tcr_file = f"tmp_for_vis_{tcr.parent.parent.id}_{tcr.id}.pdb"
            assert not pathlib.Path(tmp_tcr_file).exists()
            tmp_plip_file = "ATOM1NGLYB123_PROTEIN_UNL_Z_1.pse"
            assert not pathlib.Path(tmp_plip_file).exists()

            # test residue highlighting
            saved_session = plip_parser.create_pymol_session(
                tcr,
                save_as="./TCRpy/test/test_files/out/interactions/8gvb_test_residue_highlighted.pse",
                antigen_residues_to_highlight=[4, 6],
            )
            assert pathlib.Path(saved_session).exists()

            # test single residue highlighting
            saved_session = plip_parser.create_pymol_session(
                tcr,
                save_as="./TCRpy/test/test_files/out/interactions/8gvb_test_residue_highlighted.pse",
                antigen_residues_to_highlight=5,
            )
            assert pathlib.Path(saved_session).exists()

        except ModuleNotFoundError as e:  # tests to run if pymol is not installed
            if "pymol" not in str(e):
                raise ValueError("Only except pymol not found errors")
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                plip_parser.create_pymol_session(tcr)
                assert len(w) == 1  # check only one warning raised
                # check warning tells user to install pymol
                assert (
                    "conda install -c conda-forge -c schrodinger pymol-bundle"
                    in str(w[0].message)
                )

    def test_bound_tcr_interaction_visualisation_method(self):
        parser = TCRParser.TCRParser()
        test_file = "./TCRpy/test/test_files/8gvb.cif"
        tcr = [x for x in parser.get_tcr_structure("test_8gvb", test_file).get_TCRs()][
            0
        ]

        try:
            import pymol

            pymol_installed = True

        except ModuleNotFoundError as e:
            if "pymol" not in str(e):
                raise ValueError("Only except pymol not found errors")
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                tcr.visualise_interactions()
                assert len(w) == 1  # check only one warning raised
                # check warning tells user to install pymol
                assert (
                    "conda install -c conda-forge -c schrodinger pymol-bundle"
                    in str(w[0].message)
                )
            pymol_installed = False

        if pymol_installed:
            saved_session = tcr.visualise_interactions()
            assert pathlib.Path(saved_session).exists()
            os.remove(saved_session)  # test clean up

            # test residue highlighting
            saved_session = tcr.visualise_interactions(
                save_as="./TCRpy/test/test_files/out/interactions/8gvb_test_residue_highlighted_TCR_bound_method.pse",
                antigen_residues_to_highlight=[4, 6],
            )
            assert pathlib.Path(saved_session).exists()
