import warnings
import matplotlib.pyplot as plt
from importlib import reload
import numpy as np

from ..tcr_processing.TCRParser import TCRParser

try:
    import plip
    from plip.basic.remote import VisualizerData
    from plip.visualization.visualize import visualize_in_pymol
except ModuleNotFoundError as e:
    if "pymol" in str(e):
        warnings.warn(
            """\nPymol package not found. \nInteraction profiler initialising without visualisation capabilitites. \nTo enable pymol visualisations, install pymol with:
            \nconda install -c conda-forge -c schrodinger numpy pymol-bundle\n\n"""
        )
    elif "plip" in str(e):
        warnings.warn(
            """\n\nPLIP package not found. \nProfiling interactions will not be possible \nTo enable interaction profiling, install PLIP with:
        \npip install plip --no-deps\n\n"""
        )
except ImportError as e:
    if "pymol" in str(e):
        warnings.warn(
            f"""pymol was not imported. Interactions were not visualised. \nThis is due to an import error. Perhaps try reinstalling pymol? 
                    \nThe error trace was: {str(e)}
                    """
        )
    elif "plip" in str(e):
        warnings.warn(
            f"""\n\nPLIP was not imported. \nProfiling interactions will not be possible 
            \nThis is due to an import error. Perhaps try reinstalling plip? 
            \nThe error trace was: {str(e)}"""
        )


from . import utils as plip_utils
from .PLIPParser import PLIPParser
from .TCRpMHC_PLIP_Model_Parser import TCRpMHC_PLIP_Model_Parser


class TCRInteractionProfiler:

    def __init__(self, **kwargs):
        self.tcr_parser = TCRParser()
        self.model_parser = TCRpMHC_PLIP_Model_Parser()
        self.plip_parser = PLIPParser()

        from plip.basic import config

        config = reload(config)
        self.config = config

        if len(kwargs) > 0:
            self.set_interaction_parameters(**kwargs)

    def reset_parameters(self):
        from plip.basic import config

        config = reload(config)
        self.config = config

    def set_interaction_parameters(self, **kwargs):
        """
        Function to set global PLIP detection parameters, ie. the stcrpy API for PLIP config parameters.
        See https://github.com/pharmai/plip/blob/master/plip/plipcmd.py for how these are set in native PLIP
        See https://github.com/pharmai/plip/blob/master/plip/basic/config.py for the default values

        Default Parameters (from PLIP distribution):
            BS_DIST = 7.5  # Determines maximum distance to include binding site residues
            AROMATIC_PLANARITY = 5.0  # Determines allowed deviation from planarity in aromatic rings
            MIN_DIST = 0.5  # Minimum distance for all distance thresholds
            HYDROPH_DIST_MAX = 4.0  # Distance cutoff for detection of hydrophobic contacts
            HBOND_DIST_MAX = 4.1  # Max. distance between hydrogen bond donor and acceptor (Hubbard & Haider, 2001) + 0.6 A
            HBOND_DON_ANGLE_MIN = 100  # Min. angle at the hydrogen bond donor (Hubbard & Haider, 2001) + 10
            PISTACK_DIST_MAX = 5.5  # Max. distance for parallel or offset pistacking (McGaughey, 1998)
            PISTACK_ANG_DEV = 30  # Max. Deviation from parallel or perpendicular orientation (in degrees)
            PISTACK_OFFSET_MAX = 2.0  # Maximum offset of the two rings (corresponds to the radius of benzene + 0.5 A)
            PICATION_DIST_MAX = 6.0  # Max. distance between charged atom and aromatic ring center (Gallivan and Dougherty, 1999)
            SALTBRIDGE_DIST_MAX = 5.5  # Max. distance between centers of charge for salt bridges (Barlow and Thornton, 1983) + 1.5
            HALOGEN_DIST_MAX = 4.0  # Max. distance between oxy. and halogen (Halogen bonds in biological molecules., Auffinger)+0.5
            HALOGEN_ACC_ANGLE = 120  # Optimal acceptor angle (Halogen bonds in biological molecules., Auffinger)
            HALOGEN_DON_ANGLE = 165  # Optimal donor angle (Halogen bonds in biological molecules., Auffinger)
            HALOGEN_ANGLE_DEV = 30  # Max. deviation from optimal angle
            WATER_BRIDGE_MINDIST = 2.5  # Min. distance between water oxygen and polar atom (Jiang et al., 2005) -0.1
            WATER_BRIDGE_MAXDIST = 4.1  # Max. distance between water oxygen and polar atom (Jiang et al., 2005) +0.5
            WATER_BRIDGE_OMEGA_MIN = 71  # Min. angle between acceptor, water oxygen and donor hydrogen (Jiang et al., 2005) - 9
            WATER_BRIDGE_OMEGA_MAX = 140  # Max. angle between acceptor, water oxygen and donor hydrogen (Jiang et al., 2005)
            WATER_BRIDGE_THETA_MIN = 100  # Min. angle between water oxygen, donor hydrogen and donor atom (Jiang et al., 2005)
            METAL_DIST_MAX = 3.0  # Max. distance between metal ion and interacting atom (Harding, 2001)
            MAX_COMPOSITE_LENGTH = 200  # Filter out ligands with more than 200 fragments

        Raises:
            AttributeError: Raised if parameter being modified does not exist in config
            ValueError: Raised if value being set is not permitted.
        """

        self.reset_parameters()  # reset to ensure no leaks between configurations
        for param, value in kwargs.items():
            if not hasattr(self.config, param):
                raise AttributeError(f"PLIP self.config has no parameter {param}")

            if (
                "ANGLE" in param and not 0 < value < 180
            ):  # Check value for angle thresholds
                raise ValueError(
                    "Threshold for angles need to have values within 0 and 180."
                )
            if "DIST" in param:
                if value > 10:  # Check value for distance thresholds
                    raise ValueError(
                        "Threshold for distances must not be larger than 10 Angstrom."
                    )
                elif (
                    value > self.config.BS_DIST + 1
                ):  # Dynamically adapt the search space for binding site residues
                    self.config.BS_DIST = value + 1
            setattr(self.config, param, value)
        # Check additional conditions for interdependent thresholds
        if not self.config.HALOGEN_ACC_ANGLE > self.config.HALOGEN_ANGLE_DEV:
            raise ValueError(
                "The halogen acceptor angle has to be larger than the halogen angle deviation."
            )
        if not self.config.HALOGEN_DON_ANGLE > self.config.HALOGEN_ANGLE_DEV:
            raise ValueError(
                "The halogen donor angle has to be larger than the halogen angle deviation."
            )
        if not self.config.WATER_BRIDGE_MINDIST < self.config.WATER_BRIDGE_MAXDIST:
            raise ValueError(
                "The water bridge minimum distance has to be smaller than the water bridge maximum distance."
            )
        if not self.config.WATER_BRIDGE_OMEGA_MIN < self.config.WATER_BRIDGE_OMEGA_MAX:
            raise ValueError(
                "The water bridge omega minimum angle has to be smaller than the water bridge omega maximum angle"
            )

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
                \nconda install -c conda-forge -c schrodinger numpy pymol-bundle\n\n
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
                \nconda install -c conda-forge -c schrodinger numpy pymol-bundle\n\n
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

        plt.subplots_adjust(hspace=0.5)

        for pair in TCRA_tuples:
            heatmap_a[pair[0][1], ligand_number_mapping[int(pair[1][1])]] = (
                heatmap_a[pair[0][1], ligand_number_mapping[int(pair[1][1])]] + 1
            )

        ax_alpha.imshow(heatmap_a.T, cmap="PuRd")

        for i in plot_index:
            ax_alpha.axhline(y=i - 0.5, color="blue", linewidth=1)
            ax_alpha.axhline(y=i + 0.5, color="blue", linewidth=1)
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
        ax_beta.imshow(heatmap_b.T, cmap="PuRd")
        for i in plot_index:
            ax_beta.axhline(y=i - 0.5, color="blue", linewidth=1)
            ax_beta.axhline(y=i + 0.5, color="blue", linewidth=1)
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
