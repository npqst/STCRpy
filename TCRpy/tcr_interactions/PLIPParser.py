import typing
import pandas as pd
import warnings

import plip

# if typing.TYPE_CHECKING:
#     from ..tcr_processing.TCR import abTCR, gdTCR


from . import utils as plip_utils
from .TCRpMHC_PLIP_Model_Parser import TCRpMHC_PLIP_Model_Parser

try:
    from plip.basic.remote import VisualizerData
    from plip.visualization.visualize import visualize_in_pymol
except ModuleNotFoundError:
    warnings.warn(
        """\nPymol package not found. \nInteraction profiler initialising without visualisation capabilitites. \nTo enable pymol visualisations, install pymol with: 
        \nconda install -c conda-forge -c schrodinger pymol-bundle\n\n"""
    )


class PLIPParser:

    def parse_complex(
        self,
        complex: plip.structure.preparation.PDBComplex,
        tcr_pmhc_complex: typing.Union["abTCR", "gdTCR"] = None,
        renumbering=None,
        domain_assignment=None,
    ) -> pd.DataFrame:
        all_interactions = []
        for _, interaction_set in complex.interaction_sets.items():
            for interaction in interaction_set.all_itypes:
                try:
                    all_interactions.append(plip_utils.parse_interaction(interaction))
                except NotImplementedError as e:
                    print(e)
                    continue
        interactions_df = self._interactions_to_dataframe(all_interactions)
        if renumbering is not None and len(interactions_df) > 0:
            self._renumber_interactions(interactions_df, renumbering)
        if tcr_pmhc_complex is not None:
            self._map_amino_acids_to_ligands(interactions_df, tcr_pmhc_complex)
        if domain_assignment is not None:
            self._assign_domains_to_chains(interactions_df, domain_assignment)

        return interactions_df

    def _renumber_interactions(self, interactions_df, renumbering):
        # def imgt_number_mapping(original_idx):
        #     imgt_insertions_to_number_map = {
        #                 "A": 1,
        #                 "B": 2,
        #                 "C": 3,
        #                 "D": 4,
        #                 "E": 5,
        #                 "F": 6,
        #                 "G": 7,
        #                 "H": 8,
        #                 "I": 9,
        #             }

        #     if original_idx[-1] == ' ':
        #         return original_idx[1]
        #     else:
        #         return original_idx[1] + 0.1*imgt_insertions_to_number_map[original_idx]
        interactions_df["original_numbering"] = interactions_df.apply(
            lambda x: str(renumbering[x.protein_chain][(" ", x.protein_number, " ")][1])
            + renumbering[x.protein_chain][(" ", x.protein_number, " ")][2].strip(),
            axis=1,
        )
        return interactions_df

        for chain_id, renumber in renumbering.items():
            for plip_idx, original_idx in renumber.items():
                mask = (interactions_df.protein_chain == chain_id) & (
                    interactions_df.protein_number == plip_idx[1]
                )
                if sum(mask) > 0:
                    interactions_df[mask].loc[:, "protein_number"] = (
                        str(original_idx[1]) + str(original_idx[-1])
                    ).strip()
        return interactions_df

    def _assign_domains_to_chains(self, interactions_df, domains):
        chain_to_domain_mapping = {v: k for k, v in domains.items()}

        def assign_domain(chain_id):
            if chain_id in chain_to_domain_mapping:
                return chain_to_domain_mapping[chain_id]
            else:
                return None

        interactions_df["domain"] = interactions_df.protein_chain.apply(assign_domain)

    def _interactions_to_dataframe(self, interaction_list: list) -> pd.DataFrame:
        columns = [
            "type",
            "protein_atom",
            "protein_chain",
            "protein_residue",
            "protein_number",
            "ligand_atom",
            "distance",
            "angle",
            "plip_id",
        ]

        interactions_as_tuples = [
            interaction.to_tuple() for interaction in interaction_list
        ]
        interactions = list(zip(*interactions_as_tuples))
        if len(interactions) > 0:
            interactions_as_dict = {
                columns[i]: interactions[i] for i in range(len(columns))
            }
            return pd.DataFrame(interactions_as_dict)
        else:
            return pd.DataFrame(columns=columns)

    def _map_amino_acids_to_ligands(
        self, interactions_df: pd.DataFrame, tcr_pmhc_complex: str
    ):
        parser = TCRpMHC_PLIP_Model_Parser()
        _, _, ligand_pdb, ligand_sdf = parser.parse_tcr_pmhc_complex(
            tcr_pmhc_complex, delete_tmp_files=False, renumber=False
        )
        coordinate_mapping = parser.map_amino_acids_to_ligands(ligand_pdb, ligand_sdf)
        if len(interactions_df) > 0:
            ligand_residues = interactions_df.apply(
                lambda x: coordinate_mapping[x.ligand_atom[0][0]], axis=1
            )
            interactions_df["ligand_residue"], interactions_df["ligand_number"] = map(
                list, zip(*ligand_residues)
            )
        else:  # return empty dataframe with appropriate columns
            extended_columns = list(interactions_df.columns)
            extended_columns.extend(["ligand_residue", "ligand_number"])
            interactions_df = pd.DataFrame(columns=extended_columns)
        return interactions_df

    def _visualize_interactions(self, complex: plip.structure.preparation.PDBComplex):

        from plip.basic import config

        if not config.PYMOL:
            config.PYMOL = True
        for ligand in complex.ligands:
            complex.characterize_complex(ligand)
            visualizer_complexes = [
                VisualizerData(complex, site)
                for site in sorted(complex.interaction_sets)
                if not len(complex.interaction_sets[site].interacting_res) == 0
            ]
            try:
                visualize_in_pymol(visualizer_complexes[0])
            except NameError as e:
                warnings.warn(
                    f"""Interactions could not be visualised. Raised error {e}.
                \nTo enable pymol visualisations please install pymol in a conda environment with:
                \nconda install -c conda-forge -c schrodinger pymol-bundle\n\n
                """
                )
            return

    def create_pymol_session(
        self,
        tcr_pmhc: "TCR",
        save_as=None,
        antigen_residues_to_highlight=None,
    ):

        try:
            import pymol
            from pymol import cmd
        except ImportError as e:
            warnings.warn(
                f"""pymol could not be imported. Raised error: {str(e)}.
                \nTo enable pymol visualisations please install pymol in a conda environment with:
                \nconda install -c conda-forge -c schrodinger pymol-bundle\n\n
                """
            )
            return

        import os
        import re

        pymol.finish_launching(["pymol", "-qc"])

        plip_parser = PLIPParser()
        model_parser = TCRpMHC_PLIP_Model_Parser()

        mol = model_parser.parse_tcr_pmhc_complex(
            tcr_pmhc, renumber=True, delete_tmp_files=False
        )
        mol, prot_fil, lig_file, _, _, _ = mol
        mol.analyze()
        try:
            plip_parser.parse_complex(mol)
            plip_parser._visualize_interactions(mol)
        except (
            pymol.CmdException
        ):  # for some reason sometimes this only works the second time? Probably to do with latency in pymol loading and object selection
            plip_parser.parse_complex(mol)
            plip_parser._visualize_interactions(mol)

        pymol_session = next(
            (
                f
                for f in os.listdir(".")
                if re.match(rf"^{mol.pymol_name.upper()}.*\.pse$", f)
            ),
            None,
        )
        cmd.load(pymol_session)

        # create temporary file containing the TCR and its MHC and antigen.
        from ..tcr_processing import TCRIO

        tcrio = TCRIO.TCRIO()
        tmp_file = f"tmp_for_vis_{tcr_pmhc.parent.parent.id}_{tcr_pmhc.id}.pdb"
        tcrio.save(tcr_pmhc, save_as=tmp_file)
        cmd.load(tmp_file)

        if len(tcr_pmhc.antigen) == 1:
            antigen_chain = tcr_pmhc.antigen[0].id
            cmd.show("sticks", f"chain {antigen_chain}")
            cmd.hide("cartoon", f"chain {antigen_chain}")

            if antigen_residues_to_highlight is not None:
                if isinstance(antigen_residues_to_highlight, int):
                    antigen_residues_to_highlight = [antigen_residues_to_highlight]
                for res_nr in antigen_residues_to_highlight:
                    cmd.color(
                        "red",
                        f"chain {antigen_chain} and res {str(res_nr)} and elem C",
                    )
        else:
            if len(tcr_pmhc.antigen) == 0:
                warnings.warn(
                    f"""Could not highlight antigen, no antigen found for TCR {tcr_pmhc.parent.parent.id}_{tcr_pmhc.id}"""
                )
            else:
                warnings.warn(
                    f"""Could not highlight antigen, multiple antigen {tcr_pmhc.antigen} found for TCR {tcr_pmhc.parent.parent.id}_{tcr_pmhc.id}"""
                )

        if save_as is None:
            save_as = f"{tcr_pmhc.parent.parent.id}_{tcr_pmhc.id}_interactions.pse"

        # cmd.save(pymol_session)
        cmd.save(save_as)
        cmd.delete("all")

        # clean up pymol environment and remove temporary files
        del cmd
        os.remove(pymol_session)
        os.remove(tmp_file)

        return save_as
