"""
Created on 3rd April 2024
Nele Quast based on work by Dunbar and Leem
The TCR class.
"""

import sys
import warnings

from .Entity import Entity

try:
    from .. import tcr_interactions
except ImportError as e:
    warnings.warn(
        "TCR interaction profiling could not be imported. Check PLIP installation"
    )
    print(e)


class TCR(Entity):
    """
    TCR class. Inherits from PDB.Entity.
    This is a base class for TCR strucutres, enabling antigen and MHC association.
    abTCR and gdTCR are the instantiated subclasses of this class.
    """

    def _add_antigen(self, antigen=None):
        """
        Append associated antigen to TCR antigen field.

        Args:
            antigen (Antigen, optional): Antigen to associate with TCR. Defaults to None.
        """
        if antigen not in self.antigen:
            self.antigen.append(antigen)

    def _add_mhc(self, mhc=None):
        """
        Append associated MHC to TCR MHC field. If antigen are associted with MHC but not TCR, add them to TCR antigen.

        Args:
            mhc (MHC, optional): MHC to associate with TCR. Defaults to None.
        """
        if mhc not in self.MHC:
            self.MHC.append(mhc)
            # If there are any het antigens that are in the MHC but not in close proximity of the TCR
            # (e.g. 4x6c antigen) then add it to the TCR.
            if set(mhc.antigen) - set(self.antigen):
                self.antigen.extend(mhc.antigen)

    def get_antigen(self):
        """
        Return a list of TCR associated antigens.
        """
        return self.antigen

    def get_MHC(self):
        """
        Return a list of TCR associated MHCs.
        """
        return self.MHC

    def is_bound(self):
        """True or False if the TCR is associated with an antigen.

        Returns:
            bool: Whether TCR is associated with an antigen.
        """
        if self.get_antigen():
            return True
        else:
            return False

    def get_chains(self):
        """Returns generator of TCR chains.

        Yields:
            Chain: TCR chain
        """
        for c in self:
            yield c

    def get_residues(self):
        """Returns generator of TCR residues.

        Yields:
            Residue: TCR residue
        """
        for c in self.get_chains():
            for r in c:
                yield r

    def get_atoms(self):
        """Returns generator of TCR atoms.

        Yields:
            Atom: TCR atoms
        """
        for r in self.get_residues():
            for a in r:
                yield a

    def get_frameworks(self):
        """
        Obtain framework regions from a TCR structure object as generator.

        Yields:
            Fragment: TCR framework regions
        """
        for f in self.get_fragments():
            if "fw" in f.id:
                yield f

    def get_CDRs(self):
        """
        Obtain complementarity determining regions (CDRs) from a TCR structure object as generator.

        Yields:
            Fragment: TCR CDR regions
        """
        for f in self.get_fragments():
            if "cdr" in f.id:
                yield f

    def get_TCR_type(self):
        """Get TCR type according to variable region assignments.

        Returns:
            str: TCR type (abTCR, gdTCR, dbTCR)
        """
        if hasattr(self, "tcr_type"):
            return self.tcr_type
        elif hasattr(self, "VB") and hasattr(self, "VA"):
            self.tcr_type = "abTCR"
            return self.tcr_type
        elif hasattr(self, "VD") and hasattr(self, "VG"):
            self.tcr_type = "gdTCR"
            return self.tcr_type
        elif hasattr(self, "VB") and hasattr(self, "VD"):
            self.tcr_type = "dbTCR"
            return self.tcr_type

    def get_germline_assignments(self):
        """Retrive germline assignments for all TCR chains.
        This is a dictionary with the chain ID as key and the germline assignments as value.

        Returns:
            dict: dict with TCR chain ID as key and germline assignments as value
        """
        return {c.id: c.get_germline_assignments() for c in self.get_chains()}

    def get_MHC_allele_assignments(self):
        """
        Retrieve MHC allele assignments for all TCR associated MHCs.
        This is a list of dictionaries with the MHC ID as key and the allele assignments as value.

        Returns:
            dict: dict with MHC chain ID as key and allele assignments as value
        """
        return [
            (
                mhc.get_allele_assignments()
                if mhc.level
                != "C"  # results in identical nesting structure for MHC and MHCchain types
                else {mhc.id: mhc.get_allele_assignments()}
            )
            for mhc in self.get_MHC()
        ]

    def get_germlines_and_alleles(self):
        """Get all germline and allele assignments for TCR and MHC chains as a dictionary with the chain ID as key and the germline assignments as value.

        Returns:
            dict: Dictionary of TCR germline and MHC allele assignemnts with amino acid sequences.
        """
        from ..tcr_formats.tcr_formats import get_sequences

        germlines_and_alleles = {}

        try:
            germlines = self.get_germline_assignments()
            for tcr_domain, c in self.get_domain_assignment().items():
                germlines_and_alleles[tcr_domain] = (
                    germlines[c]["v_gene"][0][1],
                    germlines[c]["j_gene"][0][1],
                )
                germlines_and_alleles[f"{tcr_domain}_species"] = sorted(
                    tuple(
                        set(
                            (
                                germlines[c]["v_gene"][0][0],
                                germlines[c]["j_gene"][0][0],
                            )
                        )
                    )
                )
                germlines_and_alleles[f"TCR_{tcr_domain}_seq"] = get_sequences(self[c])[
                    c
                ]
            if len(self.get_MHC()) == 1:
                MHC = self.get_MHC()[0]
                alleles = self.get_MHC_allele_assignments()[0]
                germlines_and_alleles["MHC_type"] = (
                    MHC.get_MHC_type() if MHC.level != "C" else MHC.chain_type
                )
                MHC_domains = {list(d.keys())[0]: c for c, d in alleles.items()}
                for d, c in MHC_domains.items():
                    germlines_and_alleles[f"MHC_{d}"] = alleles[c][d][0][1]
                    germlines_and_alleles[f"MHC_{d}_seq"] = (
                        get_sequences(MHC[c])[c]
                        if MHC.level != "C"
                        else get_sequences(MHC)[c]
                    )
            germlines_and_alleles["antigen"] = (
                get_sequences(self.get_antigen()[0])[self.get_antigen()[0].id]
                if len(self.get_antigen()) == 1
                else None
            )
        except Exception as e:
            warnings.warn(
                f"Germline and allele retrieval failed for {self} with error {str(e)}"
            )

        return germlines_and_alleles

    def save(self, save_as=None, tcr_only: bool = False, format: str = "pdb"):
        """Save TCR object as PDB or MMCIF file.

        Args:
            save_as (str, optional): File path to save TCR to. Defaults to None.
            tcr_only (bool, optional): Whether to save TCR only or to include MHC and antigen. Defaults to False.
            format (str, optional): Whether to save as PDB or MMCIF. Defaults to "pdb".
        """
        from . import TCRIO

        tcrio = TCRIO.TCRIO()
        tcrio.save(self, save_as=save_as, tcr_only=tcr_only, format=format)

    def get_scanning_angle(self, mode="rudolph"):
        """
        Returns TCR:pMHC complex scanning (aka crossing, incident) angle of TCR to MHC.
        See paper for details.

        Args:
            mode (str, optional): Mode for calculating the scanning angle. Options "rudolph", "cys", "com". Defaults to "rudolph".

        Returns:
            float: Scanning angle of TCR to MHC in degrees
        """
        if not hasattr(self, "geometry") or self.geometry.mode != mode:
            self.calculate_docking_geometry(mode=mode)
        return self.geometry.get_scanning_angle()

    def get_pitch_angle(self, mode="cys"):
        """
        Returns TCR:pMHC complex pitch angle of TCR to MHC.
        See paper for details.

        Args:
            mode (str, optional): Mode for calculating the scanning angle. Options "rudolph", "cys", "com". Defaults to "cys".

        Returns:
            float: Pitch angle of TCR to MHC in degrees
        """
        if not hasattr(self, "geometry") or self.geometry.mode != mode:
            self.calculate_docking_geometry(mode=mode)
        return self.geometry.get_pitch_angle()

    def calculate_docking_geometry(self, mode="rudolph", as_df=False):
        """Calculate docking geometry of TCR to MHC.
        This is a wrapper function for the TCRGeom class.

        Args:
            mode (str, optional): Mode for calculating the geometry. Options "rudolph", "cys", "com". Defaults to "rudolph".
            as_df (bool, optional): Whether to return as dictionary or dataframe. Defaults to False.

        Returns:
            [dict, DataFrame]: TCR to MHC geometry.
        """
        if len(self.get_MHC()) == 0:
            warnings.warn(
                f"No MHC found for TCR {self}. Docking geometry cannot be calcuated"
            )
            return None

        try:  # import here to avoid circular imports
            from ..tcr_geometry.TCRGeom import TCRGeom
        except ImportError as e:
            warnings.warn(
                "TCR geometry calculation could not be imported. Check installation"
            )
            raise ImportError(str(e))

        self.geometry = TCRGeom(self, mode=mode)
        if as_df:
            return self.geometry.to_df()
        return self.geometry.to_dict()

    def score_docking_geometry(self, **kwargs):
        """
        Score docking geometry of TCR to MHC.
        This is a wrapper function for the TCRGeomFiltering class.
        The score is calculated as the negative log of the TCR:pMHC complex geometry feature probabilities based on the distributions fit by maximum likelihood estimation of TCR to Class I MHC strucutres from STCRDab.
        Please see the paper methods for details.

        Returns:
            float: TCR:pMHC complex score as negative log of TCR:pMHC complex geometry feature probabilities
        """
        from ..tcr_geometry.TCRGeomFiltering import DockingGeometryFilter

        geom_filter = DockingGeometryFilter()
        if not hasattr(self, "geometry"):
            self.calculate_docking_geometry(mode="com")
        return geom_filter.score_docking_geometry(
            self.geometry.get_scanning_angle(mode="com"),
            self.geometry.get_pitch_angle(mode="com"),
            self.geometry.tcr_com[-1],  # z component of TCR centre of mass
        )

    def profile_peptide_interactions(
        self, renumber: bool = True, save_to: str = None, **kwargs
    ) -> "pd.DataFrame":
        """
        Profile the interactions of the peptide to the TCR and the MHC.

        Args:
            renumber (bool, optional): Whether to renumber the interacting residues. Defaults to True.
            save_to (str, optional): Path to save intraction data to as csv. Defaults to None.

        Returns:
            pd.DataFrame: Dataframe of peptide interactions
        """
        if len(self.get_antigen()) == 0:
            warnings.warn(
                f"No peptide antigen found for TCR {self}. Peptide interactions cannot be profiled"
            )
            return None

        if "PLIPParser" not in [m.split(".")[-1] for m in sys.modules]:
            warnings.warn(
                "TCR interactions module was not imported. Check warning log and PLIP installation"
            )
            return None

        from ..tcr_interactions import TCRInteractionProfiler

        interaction_profiler = TCRInteractionProfiler.TCRInteractionProfiler(**kwargs)
        interactions = interaction_profiler.get_interactions(
            self, renumber=renumber, save_as_csv=save_to
        )
        return interactions

    def get_interaction_heatmap(self, plotting_kwargs={}, **interaction_kwargs):
        """
        Get interaction heatmap of TCR to MHC and peptide.
        Generates heatmap image.
        Plotting kwargs are passed to heatmap function.

        Args:
            plotting_kwargs (dict, optional):
                save_as: path to save heatmap image to
                interaction_type: type of interaction (eg. saltbridge, h_bond) to plot. All interactions are plotted by default.
                antigen_name: name of antigen for plot title
                mutation_index: index of antigen residues to highlight in plot
                Defaults to {
                    save_as:None,
                    interaction_type:None,
                    antigen_name:None,
                    mutation_index:None
                    }.
            interaction_kwargs: kwargs for TCRInteractionProfiler class. See TCRInteractionProfiler for details.
        """
        from ..tcr_interactions import TCRInteractionProfiler

        interaction_profiler = TCRInteractionProfiler.TCRInteractionProfiler(
            **interaction_kwargs
        )
        interaction_profiler.get_interaction_heatmap(self, **plotting_kwargs)

    def profile_TCR_interactions(self):
        raise NotImplementedError

    def profile_MHC_interactions(self):
        raise NotImplementedError

    def _create_interaction_visualiser(self):
        """Function called during TCR initialisation checks if pymol is installed and assigns a visualisation method accordingly.
        If pymol is installed, method to generate interaction visualisations is returned.
        If pymol is not installed, calling the visualisation


        Returns:
            callable: TCR bound method to visualise interactions of the TCR and MHC to the peptide.
        """
        try:
            import pymol

            def visualise_interactions(
                save_as=None, antigen_residues_to_highlight=None, **interaction_kwargs
            ):
                """Visualise peptide interactions in pymol.

                Args:
                    save_as (str, optional): path to save pymol session to. Defaults to None.
                    antigen_residues_to_highlight (list[int], optional): antigen residues to highlight red in pymol session. Defaults to None.
                    **interaction_kwargs: kwargs for TCRInteractionProfiler class. See TCRInteractionProfiler for details.
                Returns:
                    str: path to saved pymol session
                """
                from ..tcr_interactions import TCRInteractionProfiler

                interaction_profiler = TCRInteractionProfiler.TCRInteractionProfiler(
                    **interaction_kwargs
                )
                interaction_session_file = interaction_profiler.create_pymol_session(
                    self,
                    save_as=save_as,
                    antigen_residues_to_highlight=antigen_residues_to_highlight,
                )

                return interaction_session_file

            return visualise_interactions

        except ModuleNotFoundError:

            def visualise_interactions(**interaction_kwargs):
                warnings.warn(
                    f"""pymol was not imported. Interactions were not visualised.
                    \nTo enable pymol visualisations please install pymol in a conda environment with:
                    \nconda install -c conda-forge -c schrodinger numpy pymol-bundle\n\n
                    """
                )

            return visualise_interactions

        except ImportError as e:

            def visualise_interactions(import_error=e, **interaction_kwargs):
                warnings.warn(
                    f"""pymol was not imported. Interactions were not visualised. This is due to an import error. Perhaps try reinstalling pymol? 
                    \nThe error trace was: {str(import_error)}
                    """
                )

            return visualise_interactions


class abTCR(TCR):
    """
    abTCR class. Inherits from TCR.
    This is a subclass of TCR for TCRs with alpha and beta chains.
    """
    def __init__(self, c1, c2):
        """
        Initialise abTCR object. This is a subclass of TCR for TCRs with alpha and beta chains.

        Args:
            c1 (TCRchain): alpha or beta type TCR chain
            c2 (TCRchain): alpha or beta type TCR chain
        """

        if c1.chain_type == "B":
            Entity.__init__(self, c1.id + c2.id)
        else:
            Entity.__init__(self, c2.id + c1.id)

        # The TCR is a Holder class
        self.level = "H"
        self._add_domain(c1)
        self._add_domain(c2)
        self.child_list = sorted(
            self.child_list, key=lambda x: x.chain_type, reverse=True
        )  # make sure that the list goes B->A or G->D
        self.antigen = []
        self.MHC = []
        self.engineered = False
        self.scTCR = False  # This is rare but does happen

        self.visualise_interactions = self._create_interaction_visualiser()

    def __repr__(self):
        """
        String representation of the abTCR object.

        Returns:
            str: String representation of the abTCR objec
        """
        return "<TCR %s%s beta=%s; alpha=%s>" % (self.VB, self.VA, self.VB, self.VA)

    def _add_domain(self, chain):
        """
        Add a variable alpha or variable beta domain to the TCR object.
        Links the domain to the chain ID.

        Args:
            chain (TCRchain): TCR chain whose domain is being added.
        """
        if chain.chain_type == "B":
            self.VB = chain.id
        elif chain.chain_type == "A" or chain.chain_type == "D":
            self.VA = chain.id

        # Add the chain as a child of this entity.
        self.add(chain)

    def get_VB(self):
        """
        Retrieve the variable beta chain of the TCR

        Returns:
            TCRchain: VB chain
        """
        if hasattr(self, "VB"):
            return self.child_dict[self.VB]

    def get_VA(self):
        """
        Retrieve the variable alpha chain of the TCR

        Returns:
            TCRchain: VA chain
        """
        if hasattr(self, "VA"):
            return self.child_dict[self.VA]

    def get_domain_assignment(self):
        """
        Retrieve the domain assignment of the TCR as a dict with variable domain type as key and chain ID as value.

        Returns:
            dict: domain assignment from domain to chain ID, e.g. {"VA": "D", "VB": "E"}
        """
        try:
            return {"VA": self.VA, "VB": self.VB}
        except AttributeError:
            if hasattr(self, "VB"):
                return {"VB": self.VB}
            if hasattr(self, "VA"):
                return {"VA": self.VA}
        return None

    def is_engineered(self):
        """
        Flag for engineered TCRs.

        Returns:
            bool: Flag for engineered TCRs
        """
        if self.engineered:
            return True
        else:
            vb, va = self.get_VB(), self.get_VA()
            for var_domain in [vb, va]:
                if var_domain and var_domain.is_engineered():
                    self.engineered = True
                    return self.engineered

            self.engineered = False
            return False

    def get_fragments(self):
        """
        Retrieve the fragments, ie FW and CDR loops of the TCR as a generator.

        Yields:
            Fragment: fragment of TCR chain.
        """
        vb, va = self.get_VB(), self.get_VA()

        # If a variable domain exists
        for var_domain in [vb, va]:
            if var_domain:
                for frag in var_domain.get_fragments():
                    yield frag


class gdTCR(TCR):

    def __init__(self, c1, c2):

        if c1.chain_type == "D":
            Entity.__init__(self, c1.id + c2.id)
        else:
            Entity.__init__(self, c2.id + c1.id)

        # The TCR is a Holder class
        self.level = "H"
        self._add_domain(c1)
        self._add_domain(c2)
        self.child_list = sorted(
            self.child_list, key=lambda x: x.chain_type
        )  # make sure that the list goes B->A or D->G
        self.antigen = []
        self.MHC = []
        self.engineered = False
        self.scTCR = False  # This is rare but does happen

        self.visualise_interactions = self._create_interaction_visualiser()

    def __repr__(self):
        return "<TCR %s%s delta=%s; gamma=%s>" % (self.VD, self.VG, self.VD, self.VG)

    def _add_domain(self, chain):
        if chain.chain_type == "D":
            self.VD = chain.id
        elif chain.chain_type == "G":
            self.VG = chain.id

        # Add the chain as a child of this entity.
        self.add(chain)

    def get_VD(self):
        if hasattr(self, "VD"):
            return self.child_dict[self.VD]

    def get_VG(self):
        if hasattr(self, "VG"):
            return self.child_dict[self.VG]

    def get_domain_assignment(self):
        try:
            return {"VG": self.VG, "VD": self.VD}
        except AttributeError:
            if hasattr(self, "VD"):
                return {"VD": self.VD}
            if hasattr(self, "VG"):
                return {"VG": self.VG}
        return None

    def is_engineered(self):
        if self.engineered:
            return True
        else:
            vd, vg = self.get_VD(), self.get_VG()
            for var_domain in [vd, vg]:
                if var_domain and var_domain.is_engineered():
                    self.engineered = True
                    return self.engineered

            self.engineered = False
            return False

    def get_fragments(self):
        vd, vg = self.get_VD(), self.get_VG()

        # If a variable domain exists
        for var_domain in [vg, vd]:
            if var_domain:
                for frag in var_domain.get_fragments():
                    yield frag


class dbTCR(TCR):
    def __init__(self, c1, c2):
        super(TCR, self).__init__()

        if c1.chain_type == "B":
            Entity.__init__(self, c1.id + c2.id)
        else:
            Entity.__init__(self, c2.id + c1.id)

        # The TCR is a Holder class
        self.level = "H"
        self._add_domain(c1)
        self._add_domain(c2)
        self.child_list = sorted(
            self.child_list, key=lambda x: x.chain_type, reverse=False
        )  # make sure that the list goes B->D
        self.antigen = []
        self.MHC = []
        self.engineered = False
        self.scTCR = False  # This is rare but does happen

        self.visualise_interactions = self._create_interaction_visualiser()

    def __repr__(self):
        return "<TCR %s%s beta=%s; delta=%s>" % (self.VB, self.VD, self.VB, self.VD)

    def _add_domain(self, chain):
        if chain.chain_type == "B":
            self.VB = chain.id
        elif chain.chain_type == "D":
            self.VD = chain.id

        # Add the chain as a child of this entity.
        self.add(chain)

    def get_VB(self):
        if hasattr(self, "VB"):
            return self.child_dict[self.VB]

    def get_VD(self):
        if hasattr(self, "VD"):
            return self.child_dict[self.VD]

    def get_domain_assignment(self):
        try:
            return {"VD": self.VD, "VB": self.VB}
        except AttributeError:
            if hasattr(self, "VB"):
                return {"VB": self.VB}
            if hasattr(self, "VD"):
                return {"VD": self.VD}
        return None

    def is_engineered(self):
        if self.engineered:
            return True
        else:
            vb, vd = self.get_VB(), self.get_VD()
            for var_domain in [vb, vd]:
                if var_domain and var_domain.is_engineered():
                    self.engineered = True
                    return self.engineered

            self.engineered = False
            return False

    def get_fragments(self):
        vb, vd = self.get_VB(), self.get_VD()

        # If a variable domain exists
        for var_domain in [vb, vd]:
            if var_domain:
                for frag in var_domain.get_fragments():
                    yield frag
