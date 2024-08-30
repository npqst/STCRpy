import warnings
import matplotlib.pyplot as plt
import numpy as np

from ..tcr_processing.TCRParser import TCRParser

try:
    import plip
except ModuleNotFoundError:
    warnings.warn(
        """\n\nPLIP package not found. \nProfiling interactions will not be possible \nTo enable interaction profiling, install PLIP with:
        \npip install plip --no-deps\n\n"""
    )


from . import utils as plip_utils
from .PLIPParser import PLIPParser
from .TCRpMHC_PLIP_Model_Parser import TCRpMHC_PLIP_Model_Parser


try:
    from plip.basic.remote import VisualizerData
    from plip.visualization.visualize import visualize_in_pymol
except ModuleNotFoundError as e:
    if "pymol" in str(e):
        warnings.warn(
            """\nPymol package not found. \nInteraction profiler initialising without visualisation capabilitites. \nTo enable pymol visualisations, install pymol with:
            \nconda install -c conda-forge -c schrodinger numpy==1.26.0 pymol-bundle\n\n"""
        )
    elif "plip" in str(e):
        warnings.warn(
            """\n\nPLIP package not found. \nProfiling interactions will not be possible \nTo enable interaction profiling, install PLIP with:
        \npip install plip --no-deps\n\n"""
        )


class TCRInteractionProfiler:
    def __init__(self):
        self.tcr_parser = TCRParser()
        self.model_parser = TCRpMHC_PLIP_Model_Parser()
        self.plip_parser = PLIPParser()

    def _visualize_interactions(self, complex: "plip.structure.preparation.PDBComplex"):

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
                \nconda install -c conda-forge -c schrodinger numpy==1.26.0 pymol-bundle\n\n
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
                \nconda install -c conda-forge -c schrodinger numpy==1.26.0 pymol-bundle\n\n
                """
            )
            return

        import os
        import re

        pymol.finish_launching(["pymol", "-qc"])

        mol = self.model_parser.parse_tcr_pmhc_complex(
            tcr_pmhc, renumber=True, delete_tmp_files=True
        )
        mol, _, _ = mol
        mol.analyze()
        try:
            self.plip_parser.parse_complex(mol)
            self._visualize_interactions(mol)
        except (
            pymol.CmdException
        ):  # for some reason sometimes this only works the second time? Probably to do with latency in pymol loading and object selection
            self.plip_parser.parse_complex(mol)
            self._visualize_interactions(mol)

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

    def get_interactions(self, tcr, renumber=True, save_as_csv=None):
        mol = self.model_parser.parse_tcr_pmhc_complex(tcr, renumber=renumber)
        if renumber:
            mol, renumbering, domains = mol
        else:
            renumbering = None
            domains = None
        mol.analyze()

        interactions_df = self.plip_parser.parse_complex(
            mol, tcr, renumbering=renumbering, domain_assignment=domains
        )

        if save_as_csv is not None:
            interactions_df.to_csv(save_as_csv)

        return interactions_df

    def get_interaction_heatmap(self, tcr, renumber=True, **plotting_kwargs):
        interactions_df = self.get_interactions(tcr, renumber=renumber)

        heatmaps = self._interaction_heatmap(
            interactions_df,
            tcr_name=f"{tcr.parent.parent.id}_{tcr.id}",
            peptide_length=len(tcr.antigen[0]),
            **plotting_kwargs,
        )
        return heatmaps

    @staticmethod
    def _interaction_heatmap(
        interactions_df,
        tcr_name=None,
        peptide_length=10,
        save_as=None,
        interaction_type=None,
        antigen_name=None,
        mutation_index=None,
    ):

        if interaction_type is not None:
            df = interactions_df[interactions_df.type == interaction_type]
        else:
            df = interactions_df

        if antigen_name is None:
            antigen_name = "peptide"

        TCRA_interactions = df[df.domain.apply(lambda x: x in ["VA", "VD"])]
        TCRB_interactions = df[df.domain == "VB"]
        TCRA_tuples = TCRA_interactions.apply(
            lambda x: (
                (x["protein_residue"], x["protein_number"]),
                (x["ligand_residue"], x["ligand_number"]),
            ),
            axis=1,
        )
        TCRB_tuples = TCRB_interactions.apply(
            lambda x: (
                (x["protein_residue"], x["protein_number"]),
                (x["ligand_residue"], x["ligand_number"]),
            ),
            axis=1,
        )

        heatmap_a = np.zeros((126, peptide_length))
        heatmap_b = np.zeros((126, peptide_length))

        # check peptide numbering
        offset = max(set(interactions_df.ligand_number)) + 1 - peptide_length
        ligand_number_mapping = {x + int(offset): x for x in range(peptide_length)}

        if "original_numbering" in interactions_df.columns:
            tcr_a_mapping = list(
                zip(
                    *set(
                        [
                            (
                                x.protein_number,
                                f"{x.original_numbering}-{x.protein_residue}",
                            )
                            for _, x in interactions_df.iterrows()
                            if x.domain in ["VA", "VD"]
                        ]
                    )
                )
            )
            tcr_b_mapping = list(
                zip(
                    *set(
                        [
                            (
                                x.protein_number,
                                f"{x.original_numbering}-{x.protein_residue}",
                            )
                            for _, x in interactions_df.iterrows()
                            if x.domain == "VB"
                        ]
                    )
                )
            )
            peptide_mapping = list(
                zip(
                    *set(
                        [
                            (
                                x.ligand_number - offset,
                                f"{x.ligand_number}-{x.ligand_residue}",
                            )
                            for _, x in interactions_df.iterrows()
                        ]
                    )
                )
            )
            peptide_mapping_dict = dict(zip(*reversed(peptide_mapping)))

        if mutation_index is not None:
            if isinstance(mutation_index, str):
                mutation_index = [mutation_index]
            try:
                plot_index = [peptide_mapping_dict[m_idx] for m_idx in mutation_index]
            except KeyError:
                plot_index = []
                warnings.warn(
                    f"Mutation index could not be resolved. Peptide residues are: {list(peptide_mapping_dict.keys())}"
                )

        else:
            plot_index = []

        if interaction_type is None:
            interaction_type = "all"

        fig, (ax_alpha, ax_beta) = plt.subplots(2, 1, figsize=(18, 4))

        for pair in TCRA_tuples:
            heatmap_a[pair[0][1], ligand_number_mapping[int(pair[1][1])]] = (
                heatmap_a[pair[0][1], ligand_number_mapping[int(pair[1][1])]] + 1
            )

        ax_alpha.imshow(heatmap_a.T)
        for i in plot_index:
            ax_alpha.axhline(y=i - 0.5, color="red", linewidth=2)
            ax_alpha.axhline(y=i + 0.5, color="red", linewidth=2)
        ax_alpha.set_title(
            f"{tcr_name} TCR alpha chain to {antigen_name}; {interaction_type} interactions"
        )
        if len(tcr_a_mapping) > 0:
            ax_alpha.set_xticks(tcr_a_mapping[0], tcr_a_mapping[1], rotation=90)
            ax_alpha.set_yticks(peptide_mapping[0], peptide_mapping[1])
        else:
            ax_alpha.set_xticks([], [], rotation=90)
            ax_alpha.set_yticks([], [])

        for pair in TCRB_tuples:
            heatmap_b[pair[0][1], ligand_number_mapping[int(pair[1][1])]] = (
                heatmap_b[pair[0][1], ligand_number_mapping[int(pair[1][1])]] + 1
            )
        ax_beta.imshow(heatmap_b.T)
        for i in plot_index:
            ax_beta.axhline(y=i - 0.5, color="red", linewidth=2)
            ax_beta.axhline(y=i + 0.5, color="red", linewidth=2)
        ax_beta.set_title(
            f"{tcr_name} TCR beta chain to {antigen_name}; {interaction_type} interactions"
        )
        if len(tcr_b_mapping) > 0:
            ax_beta.set_xticks(tcr_b_mapping[0], tcr_b_mapping[1], rotation=90)
            ax_beta.set_yticks(peptide_mapping[0], peptide_mapping[1])
        else:
            ax_beta.set_xticks([], [], rotation=90)
            ax_beta.set_yticks([], [])

        if save_as is not None:
            fig.savefig(save_as, bbox_inches="tight", dpi=200)

        return {"alpha": heatmap_a, "beta": heatmap_b}
