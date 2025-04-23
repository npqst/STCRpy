import typing
import warnings
import os
import pandas as pd


from . import utils as plip_utils
from .TCRpMHC_PLIP_Model_Parser import TCRpMHC_PLIP_Model_Parser

class PLIPParser:
    """This class is used to parse the interactions of a TCR-pMHC complex using PLIP."""
    def parse_complex(
        self,
        complex: "plip.structure.preparation.PDBComplex",
        tcr_pmhc_complex: typing.Union["abTCR", "gdTCR"] = None,
        renumbering=None,
        domain_assignment=None,
    ) -> pd.DataFrame:
        """
        Parses PLIP profiled interactions and maps them back onto a syctpy TCR object.

        Args:
            complex (plip.structure.preparation.PDBComplex):
            tcr_pmhc_complex (typing.Union[&quot;abTCR&quot;, &quot;gdTCR&quot;], optional): _description_. Defaults to None.
            renumbering (_type_, optional): _description_. Defaults to None.
            domain_assignment (_type_, optional): _description_. Defaults to None.

        Returns:
            pd.DataFrame: _description_
        """
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
        _, tcr_mhc_pdb, ligand_pdb, ligand_sdf = parser.parse_tcr_pmhc_complex(
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

        # delete temp files needed for renumbering
        os.remove(tcr_mhc_pdb)
        os.remove(ligand_pdb)
        os.remove(ligand_sdf)

        return interactions_df
