import unittest
import pathlib
import warnings
import os

try:
    from plip.structure.preparation import PDBComplex
except ModuleNotFoundError:
    pass

from stcrpy.tcr_processing import TCRParser
from stcrpy.tcr_interactions.TCRpMHC_PLIP_Model_Parser import (
    TCRpMHC_PLIP_Model_Parser,
)
from stcrpy.tcr_interactions.PLIPParser import PLIPParser
from stcrpy.tcr_interactions.TCRInteractionProfiler import TCRInteractionProfiler


class TestTCRInteractions(unittest.TestCase):

    def test_tcrpmhc_plip_model_parser(self):

        parser = TCRParser.TCRParser()

        model_parser = TCRpMHC_PLIP_Model_Parser()

        test_file = "./test_files/8gvb.cif"
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
        test_file = "./test_files/8gvb.cif"
        tcr = [x for x in parser.get_tcr_structure("tmp", test_file).get_TCRs()][0]
        mol, renumbering, domains = model_parser.parse_tcr_pmhc_complex(tcr)
        mol.analyze()

        plip_parser = PLIPParser()
        interactions = plip_parser.parse_complex(mol, tcr, renumbering, domains)

        assert len(interactions) == 27
        assert len(interactions[interactions.type == "hbond"]) == 11
        assert len(interactions[interactions.type == "hydrophobic"]) == 12
        assert len(interactions[interactions.type == "pistack"]) == 1
        assert len(interactions[interactions.type == "saltbridge"]) == 3

        assert len(interactions[interactions.domain == "VB"]) == 1
        assert interactions[interactions.domain == "VB"].protein_residue.item() == "ASP"
        assert interactions[interactions.domain == "VB"].protein_number.item() == 110

    def test_TCR_interaction_profiler(self):
        parser = TCRParser.TCRParser()
        test_file = "./test_files/8gvb.cif"
        tcr = [x for x in parser.get_tcr_structure("tmp", test_file).get_TCRs()][0]

        interaction_profiler = TCRInteractionProfiler()
        interactions = interaction_profiler.get_interactions(tcr, renumber=True)

        assert len(interactions) == 27
        assert len(interactions[interactions.type == "hbond"]) == 11
        assert len(interactions[interactions.type == "hydrophobic"]) == 12
        assert len(interactions[interactions.type == "pistack"]) == 1
        assert len(interactions[interactions.type == "saltbridge"]) == 3

        assert len(interactions[interactions.domain == "VB"]) == 1
        assert interactions[interactions.domain == "VB"].protein_residue.item() == "ASP"
        assert interactions[interactions.domain == "VB"].protein_number.item() == 110

        interactions = interaction_profiler.get_interactions(tcr, renumber=False)
        assert len(interactions) == 27
        assert len(interactions[interactions.type == "hbond"]) == 11
        assert len(interactions[interactions.type == "hydrophobic"]) == 12
        assert len(interactions[interactions.type == "pistack"]) == 1
        assert len(interactions[interactions.type == "saltbridge"]) == 3

        csv_path = "./test_files/out/interactions/test_8gvb_interactions.csv"
        if pathlib.Path(csv_path).exists():
            os.remove(csv_path)
        interactions = interaction_profiler.get_interactions(
            tcr,
            save_as_csv=csv_path,
        )
        assert pathlib.Path(csv_path).exists()

    def test_TCR_plip_methods(self):
        parser = TCRParser.TCRParser()
        test_file = "./test_files/8gvb.cif"
        tcr = [x for x in parser.get_tcr_structure("tmp", test_file).get_TCRs()][0]

        interactions = tcr.profile_peptide_interactions()

        assert len(interactions) == 27
        assert len(interactions[interactions.type == "hbond"]) == 11
        assert len(interactions[interactions.type == "hydrophobic"]) == 12
        assert len(interactions[interactions.type == "pistack"]) == 1
        assert len(interactions[interactions.type == "saltbridge"]) == 3

        assert len(interactions[interactions.domain == "VB"]) == 1
        assert interactions[interactions.domain == "VB"].protein_residue.item() == "ASP"
        assert interactions[interactions.domain == "VB"].protein_number.item() == 110

    def test_pymol_visualisation(self):
        parser = TCRParser.TCRParser()
        model_parser = TCRpMHC_PLIP_Model_Parser()
        test_file = "./test_files/8gvb.cif"
        tcr = [x for x in parser.get_tcr_structure("tmp", test_file).get_TCRs()][0]
        mol, renumbering, domains = model_parser.parse_tcr_pmhc_complex(tcr)
        # mol.analyze()

        # plip_parser = PLIPParser()

        interaction_profiler = TCRInteractionProfiler()

        # test if plip visualisations are generated

        pymol_plip_session_name = "./ATOM1NGLYB123_PROTEIN_UNL_Z_1.pse"
        if os.path.exists(pymol_plip_session_name):
            os.remove(
                pymol_plip_session_name
            )  # remove so test can check if file is generated

        try:  # test that should complete if pymol is installed
            import pymol

            try:
                interaction_profiler._visualize_interactions(mol)
            except (
                pymol.CmdException
            ):  # sometimes function needs to run twice? Probably due to pymol loading and object selection latency
                interaction_profiler._visualize_interactions(mol)
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
                interaction_profiler._visualize_interactions(mol)
                assert len(w) == 1  # check only one warning raised
                # check warning tells user to install pymol
                assert (
                    "conda install -c conda-forge -c schrodinger numpy pymol-bundle"
                    in str(w[0].message)
                )

    def test_create_pymol_session(self):
        parser = TCRParser.TCRParser()
        test_file = "./test_files/8gvb.cif"
        tcr = [x for x in parser.get_tcr_structure("test_8gvb", test_file).get_TCRs()][
            0
        ]

        interaction_profiler = TCRInteractionProfiler()

        try:  # tests to run if pymol is installed
            import pymol

            # test automated file name generation
            saved_session = interaction_profiler.create_pymol_session(tcr)
            assert saved_session == f"{tcr.parent.parent.id}_{tcr.id}_interactions.pse"
            assert pathlib.Path(saved_session).exists()
            os.remove(saved_session)  # clean up after test

            # test saving to specified file
            session_file = "./test_files/out/interactions/8gvb_test.pse"
            saved_session = interaction_profiler.create_pymol_session(
                tcr, save_as=session_file
            )
            assert pathlib.Path(saved_session).exists()
            assert session_file == saved_session

            # check tmp file clean up works
            tmp_tcr_file = f"tmp_for_vis_{tcr.parent.parent.id}_{tcr.id}.pdb"
            assert not pathlib.Path(tmp_tcr_file).exists()
            tmp_plip_file = "ATOM1NGLYB123_PROTEIN_UNL_Z_1.pse"
            assert not pathlib.Path(tmp_plip_file).exists()

            # test residue highlighting
            saved_session = interaction_profiler.create_pymol_session(
                tcr,
                save_as="./test_files/out/interactions/8gvb_test_residue_highlighted.pse",
                antigen_residues_to_highlight=[4, 6],
            )
            assert pathlib.Path(saved_session).exists()

            # test single residue highlighting
            saved_session = interaction_profiler.create_pymol_session(
                tcr,
                save_as="./test_files/out/interactions/8gvb_test_residue_highlighted.pse",
                antigen_residues_to_highlight=5,
            )
            assert pathlib.Path(saved_session).exists()

        except ModuleNotFoundError as e:  # tests to run if pymol is not installed
            if "pymol" not in str(e):
                raise ValueError("Only except pymol not found errors")
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                interaction_profiler.create_pymol_session(tcr)
                assert len(w) == 1  # check only one warning raised
                # check warning tells user to install pymol
                assert (
                    "conda install -c conda-forge -c schrodinger numpy pymol-bundle"
                    in str(w[0].message)
                )

    def test_bound_tcr_interaction_visualisation_method(self):
        parser = TCRParser.TCRParser()
        test_file = "./test_files/8gvb.cif"
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
                    "conda install -c conda-forge -c schrodinger numpy pymol-bundle"
                    in str(w[0].message)
                )
            pymol_installed = False

        if pymol_installed:
            saved_session = tcr.visualise_interactions()
            assert pathlib.Path(saved_session).exists()
            os.remove(saved_session)  # test clean up

            # test residue highlighting
            saved_session = tcr.visualise_interactions(
                save_as="./test_files/out/interactions/8gvb_test_residue_highlighted_TCR_bound_method.pse",
                antigen_residues_to_highlight=[4, 6],
            )
            assert pathlib.Path(saved_session).exists()

    # def test_interaction_heatmap(self):           ## matplotlib pyplot kills vscode unit test suite
    #     parser = TCRParser.TCRParser()
    #     test_file = "./test_files/8gvb.cif"
    #     tcr = [x for x in parser.get_tcr_structure("test_8gvb", test_file).get_TCRs()][
    #         0
    #     ]

    #     interaction_profiler = TCRInteractionProfiler()
    #     heatmaps = interaction_profiler.get_interaction_heatmap(
    #         tcr, save_as="./examples/example_8gvb_interaction_heatmap.png"
    #     )

    def test_set_interaction_parameters(self):

        from plip.basic import config

        interaction_profiler = TCRInteractionProfiler()

        assert (
            interaction_profiler.config.HBOND_DON_ANGLE_MIN
            == 100.0
            == config.HBOND_DON_ANGLE_MIN
        )

        assert interaction_profiler.config.BS_DIST == 7.5 == config.BS_DIST

        interaction_profiler.set_interaction_parameters(
            HBOND_DON_ANGLE_MIN=25, BS_DIST=8.5
        )

        assert (
            interaction_profiler.config.HBOND_DON_ANGLE_MIN
            == 25.0
            == config.HBOND_DON_ANGLE_MIN
        )

        assert interaction_profiler.config.BS_DIST == 8.5 == config.BS_DIST

        interaction_profiler.reset_parameters()

        assert (
            interaction_profiler.config.HBOND_DON_ANGLE_MIN
            == 100.0
            == config.HBOND_DON_ANGLE_MIN
        )

        assert interaction_profiler.config.BS_DIST == 7.5 == config.BS_DIST

        # same test but set paramters from initialisation
        interaction_profiler = TCRInteractionProfiler(
            HBOND_DON_ANGLE_MIN=25, BS_DIST=8.5
        )

        assert (
            interaction_profiler.config.HBOND_DON_ANGLE_MIN
            == 25.0
            == config.HBOND_DON_ANGLE_MIN
        )

        assert interaction_profiler.config.BS_DIST == 8.5 == config.BS_DIST

        interaction_profiler.reset_parameters()

        assert (
            interaction_profiler.config.HBOND_DON_ANGLE_MIN
            == 100.0
            == config.HBOND_DON_ANGLE_MIN
        )

        assert interaction_profiler.config.BS_DIST == 7.5 == config.BS_DIST

    def test_setting_and_using_alternative_interaction_parameters(self):
        import stcrpy

        tcr = stcrpy.load_TCRs("test_files/8gvb.cif")[0]

        interaction_profiler = TCRInteractionProfiler(
            BS_DIST=10.0,
            HYDROPH_DIST_MAX=6.0,
            HBOND_DIST_MAX=5.0,
            HBOND_DON_ANGLE_MIN=1.0,
        )  # these parameters are more permissive -> leads to larger dataframe of interactions
        alt_interactions_df = interaction_profiler.get_interactions(tcr)

        default_interaction_profiler = TCRInteractionProfiler()

        default_interactions_df = default_interaction_profiler.get_interactions(tcr)

        interaction_profiler.reset_parameters()
        reset_interactions_df = interaction_profiler.get_interactions(tcr)

        assert len(alt_interactions_df) > len(default_interactions_df)
        assert len(default_interactions_df) == len(reset_interactions_df)

    def test_bound_method_alternative_interaction_parameters(self):
        import stcrpy

        tcr = stcrpy.load_TCRs("test_files/8gvb.cif")[0]

        alt_interactions_df = tcr.profile_peptide_interactions(
            BS_DIST=10.0,
            HYDROPH_DIST_MAX=6.0,
            HBOND_DIST_MAX=5.0,
            HBOND_DON_ANGLE_MIN=1.0,
        )

        default_interactions_df = tcr.profile_peptide_interactions()
        assert len(alt_interactions_df) > len(default_interactions_df)

    def test_unconventional_peptide_profiling(self):
        import stcrpy

        tcr = stcrpy.fetch_TCR("6u3n")
        interactions = tcr.profile_peptide_interactions()
        assert len(interactions) == 15

        tcr1, tcr2 = stcrpy.fetch_TCR("4pjf")
        interactions = tcr1.profile_peptide_interactions()
        assert len(interactions) == 11
        interactions = tcr2.profile_peptide_interactions()
        assert len(interactions) == 10

        tcr1, tcr2 = stcrpy.fetch_TCR("5d7i")
        interactions = tcr1.profile_peptide_interactions()
        assert len(interactions) == 10
        interactions = tcr2.profile_peptide_interactions()
        assert len(interactions) == 10

        tcr = stcrpy.fetch_TCR("3arb")
        interactions = tcr.profile_peptide_interactions()
        assert len(interactions) == 19
