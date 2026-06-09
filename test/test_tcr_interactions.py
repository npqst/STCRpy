import unittest
import pathlib
import warnings
import os

try:
    from plip.structure.preparation import PDBComplex
except ModuleNotFoundError:
    pass

import stcrpy

try:
    import plip

    HAS_PLIP = True
    from stcrpy.tcr_interactions.TCRpMHC_PLIP_Model_Parser import (
        TCRpMHC_PLIP_Model_Parser,
    )
    from stcrpy.tcr_interactions.PLIPParser import PLIPParser
    from stcrpy.tcr_interactions.TCRInteractionProfiler import TCRInteractionProfiler

except ImportError:
    HAS_PLIP = False

STRUCTURES = "./test_files/structures"


@unittest.skipUnless(HAS_PLIP, "plip is not installed")
class TestTCRInteractions(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tcr = stcrpy.load_TCRs(f"{STRUCTURES}/8gvb.cif")[0]

    def test_tcrpmhc_plip_model_parser(self):
        model_parser = TCRpMHC_PLIP_Model_Parser()
        mol, _, domains = model_parser.parse_tcr_pmhc_complex(self.tcr)
        assert isinstance(mol, PDBComplex)
        assert domains == {"VA": "A", "VB": "B"}
        mol.analyze()

        mol = model_parser.parse_tcr_pmhc_complex(self.tcr, renumber=False)
        assert isinstance(mol, PDBComplex)
        mol.analyze()

    def test_plip_parser(self):
        model_parser = TCRpMHC_PLIP_Model_Parser()
        mol, renumbering, domains = model_parser.parse_tcr_pmhc_complex(self.tcr)
        mol.analyze()

        plip_parser = PLIPParser()
        interactions = plip_parser.parse_complex(mol, self.tcr, renumbering, domains)

        assert len(interactions) == 27
        assert len(interactions[interactions.type == "hbond"]) == 11
        assert len(interactions[interactions.type == "hydrophobic"]) == 12
        assert len(interactions[interactions.type == "pistack"]) == 1
        assert len(interactions[interactions.type == "saltbridge"]) == 3

        assert len(interactions[interactions.domain == "VB"]) == 1
        assert interactions[interactions.domain == "VB"].protein_residue.item() == "ASP"
        assert interactions[interactions.domain == "VB"].protein_number.item() == 110

    def test_TCR_interaction_profiler(self):
        interaction_profiler = TCRInteractionProfiler()
        interactions = interaction_profiler.get_interactions(self.tcr, renumber=True)

        assert len(interactions) == 27
        assert len(interactions[interactions.type == "hbond"]) == 11
        assert len(interactions[interactions.type == "hydrophobic"]) == 12
        assert len(interactions[interactions.type == "pistack"]) == 1
        assert len(interactions[interactions.type == "saltbridge"]) == 3

        assert len(interactions[interactions.domain == "VB"]) == 1
        assert interactions[interactions.domain == "VB"].protein_residue.item() == "ASP"
        assert interactions[interactions.domain == "VB"].protein_number.item() == 110

        interactions = interaction_profiler.get_interactions(self.tcr, renumber=False)
        assert len(interactions) == 27
        assert len(interactions[interactions.type == "hbond"]) == 11
        assert len(interactions[interactions.type == "hydrophobic"]) == 12
        assert len(interactions[interactions.type == "pistack"]) == 1
        assert len(interactions[interactions.type == "saltbridge"]) == 3

        csv_path = "./test_files/out/interactions/test_8gvb_interactions.csv"
        if pathlib.Path(csv_path).exists():
            os.remove(csv_path)
        interactions = interaction_profiler.get_interactions(
            self.tcr,
            save_as_csv=csv_path,
        )
        assert pathlib.Path(csv_path).exists()

    def test_TCR_plip_methods(self):
        interactions = self.tcr.profile_peptide_interactions()

        assert len(interactions) == 27
        assert len(interactions[interactions.type == "hbond"]) == 11
        assert len(interactions[interactions.type == "hydrophobic"]) == 12
        assert len(interactions[interactions.type == "pistack"]) == 1
        assert len(interactions[interactions.type == "saltbridge"]) == 3

        assert len(interactions[interactions.domain == "VB"]) == 1
        assert interactions[interactions.domain == "VB"].protein_residue.item() == "ASP"
        assert interactions[interactions.domain == "VB"].protein_number.item() == 110

    def test_pymol_visualisation(self):
        model_parser = TCRpMHC_PLIP_Model_Parser()
        mol, _, _ = model_parser.parse_tcr_pmhc_complex(self.tcr)
        interaction_profiler = TCRInteractionProfiler()

        pymol_plip_session_name = "./ATOM1NGLYB323_PROTEIN_UNL_Z_1.pse"
        if os.path.exists(pymol_plip_session_name):
            os.remove(pymol_plip_session_name)

        try:
            import pymol

            try:
                interaction_profiler._visualize_interactions(mol)
            except pymol.CmdException:
                interaction_profiler._visualize_interactions(mol)

            path = pathlib.Path(pymol_plip_session_name)
            assert path.is_file()
            os.remove(pymol_plip_session_name)
        except ModuleNotFoundError as e:
            if "pymol" not in str(e):
                raise ValueError("Only except pymol not found errors")
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                interaction_profiler._visualize_interactions(mol)
                assert len(w) == 1
                assert (
                    "conda install -c conda-forge -c schrodinger numpy pymol-bundle"
                    in str(w[0].message)
                )

    def test_create_pymol_session(self):
        interaction_profiler = TCRInteractionProfiler()

        try:
            import pymol

            saved_session = interaction_profiler.create_pymol_session(self.tcr)
            assert saved_session == f"{self.tcr.parent.parent.id}_{self.tcr.id}_interactions.pse"
            assert pathlib.Path(saved_session).exists()
            os.remove(saved_session)

            session_file = "./test_files/out/interactions/8gvb_test.pse"
            saved_session = interaction_profiler.create_pymol_session(
                self.tcr, save_as=session_file
            )
            assert pathlib.Path(saved_session).exists()
            assert session_file == saved_session

            tmp_tcr_file = f"tmp_for_vis_{self.tcr.parent.parent.id}_{self.tcr.id}.pdb"
            assert not pathlib.Path(tmp_tcr_file).exists()
            tmp_plip_file = "ATOM1NGLYB123_PROTEIN_UNL_Z_1.pse"
            assert not pathlib.Path(tmp_plip_file).exists()

            saved_session = interaction_profiler.create_pymol_session(
                self.tcr,
                save_as="./test_files/out/interactions/8gvb_test_residue_highlighted.pse",
                antigen_residues_to_highlight=[4, 6],
            )
            assert pathlib.Path(saved_session).exists()

            saved_session = interaction_profiler.create_pymol_session(
                self.tcr,
                save_as="./test_files/out/interactions/8gvb_test_residue_highlighted.pse",
                antigen_residues_to_highlight=5,
            )
            assert pathlib.Path(saved_session).exists()

        except ModuleNotFoundError as e:
            if "pymol" not in str(e):
                raise ValueError("Only except pymol not found errors")
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                interaction_profiler.create_pymol_session(self.tcr)
                assert len(w) == 1
                assert (
                    "conda install -c conda-forge -c schrodinger numpy pymol-bundle"
                    in str(w[0].message)
                )

    def test_bound_tcr_interaction_visualisation_method(self):
        try:
            import pymol

            pymol_installed = True
        except ModuleNotFoundError as e:
            if "pymol" not in str(e):
                raise ValueError("Only except pymol not found errors")
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                self.tcr.visualise_interactions()
                assert len(w) == 1
                assert (
                    "conda install -c conda-forge -c schrodinger numpy pymol-bundle"
                    in str(w[0].message)
                )
            pymol_installed = False

        if pymol_installed:
            saved_session = self.tcr.visualise_interactions()
            assert pathlib.Path(saved_session).exists()
            os.remove(saved_session)

            saved_session = self.tcr.visualise_interactions(
                save_as="./test_files/out/interactions/8gvb_test_residue_highlighted_TCR_bound_method.pse",
                antigen_residues_to_highlight=[4, 6],
            )
            assert pathlib.Path(saved_session).exists()

    def test_set_interaction_parameters(self):
        from plip.basic import config

        interaction_profiler = TCRInteractionProfiler()

        assert interaction_profiler.config.HBOND_DON_ANGLE_MIN == 100.0 == config.HBOND_DON_ANGLE_MIN
        assert interaction_profiler.config.BS_DIST == 7.5 == config.BS_DIST

        interaction_profiler.set_interaction_parameters(HBOND_DON_ANGLE_MIN=25, BS_DIST=8.5)

        assert interaction_profiler.config.HBOND_DON_ANGLE_MIN == 25.0 == config.HBOND_DON_ANGLE_MIN
        assert interaction_profiler.config.BS_DIST == 8.5 == config.BS_DIST

        interaction_profiler.reset_parameters()

        assert interaction_profiler.config.HBOND_DON_ANGLE_MIN == 100.0 == config.HBOND_DON_ANGLE_MIN
        assert interaction_profiler.config.BS_DIST == 7.5 == config.BS_DIST

        interaction_profiler = TCRInteractionProfiler(HBOND_DON_ANGLE_MIN=25, BS_DIST=8.5)

        assert interaction_profiler.config.HBOND_DON_ANGLE_MIN == 25.0 == config.HBOND_DON_ANGLE_MIN
        assert interaction_profiler.config.BS_DIST == 8.5 == config.BS_DIST

        interaction_profiler.reset_parameters()

        assert interaction_profiler.config.HBOND_DON_ANGLE_MIN == 100.0 == config.HBOND_DON_ANGLE_MIN
        assert interaction_profiler.config.BS_DIST == 7.5 == config.BS_DIST

    def test_setting_and_using_alternative_interaction_parameters(self):
        interaction_profiler = TCRInteractionProfiler(
            BS_DIST=10.0,
            HYDROPH_DIST_MAX=6.0,
            HBOND_DIST_MAX=5.0,
            HBOND_DON_ANGLE_MIN=1.0,
        )
        alt_interactions_df = interaction_profiler.get_interactions(self.tcr)

        default_interaction_profiler = TCRInteractionProfiler()
        default_interactions_df = default_interaction_profiler.get_interactions(self.tcr)

        interaction_profiler.reset_parameters()
        reset_interactions_df = interaction_profiler.get_interactions(self.tcr)

        assert len(alt_interactions_df) > len(default_interactions_df)
        assert len(default_interactions_df) == len(reset_interactions_df)

    def test_bound_method_alternative_interaction_parameters(self):
        alt_interactions_df = self.tcr.profile_peptide_interactions(
            BS_DIST=10.0,
            HYDROPH_DIST_MAX=6.0,
            HBOND_DIST_MAX=5.0,
            HBOND_DON_ANGLE_MIN=1.0,
        )
        default_interactions_df = self.tcr.profile_peptide_interactions()
        assert len(alt_interactions_df) > len(default_interactions_df)

    def test_unconventional_peptide_profiling(self):
        tcr = stcrpy.load_TCRs(f"{STRUCTURES}/6u3n.cif")[0]
        interactions = tcr.profile_peptide_interactions()
        assert len(interactions) == 15

        tcr1, tcr2 = stcrpy.load_TCRs(f"{STRUCTURES}/4pjf.cif")
        interactions = tcr1.profile_peptide_interactions()
        assert len(interactions) == 11
        interactions = tcr2.profile_peptide_interactions()
        assert len(interactions) == 10

        tcr1, tcr2 = stcrpy.load_TCRs(f"{STRUCTURES}/5d7i.cif")
        interactions = tcr1.profile_peptide_interactions()
        assert len(interactions) == 10
        interactions = tcr2.profile_peptide_interactions()
        assert len(interactions) == 10

        tcr = stcrpy.load_TCRs(f"{STRUCTURES}/3arb.cif")[0]
        interactions = tcr.profile_peptide_interactions()
        assert len(interactions) == 19
