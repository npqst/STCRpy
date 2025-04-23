import os
import Bio
from typing import Union
import warnings
import numpy as np

from ..tcr_processing.TCRParser import TCRParser
from ..tcr_processing.TCRIO import TCRIO
from ..tcr_processing import abTCR, MHCchain


# Some of this code is adapted and refactored from https://github.com/EsamTolba/TCR-CoM/
# Please see the TCRCoM_LICENSE for the license that applies to those code sections.

class TCRCoM:
    def __init__(self):
        """Abstract class for calculating TCR centre of mass after aligning TCR:pMHC complex to reference MHC structure."""
        self.set_reffile()

        tcr_parser = TCRParser()
        self.ref_model = list(
            tcr_parser.get_tcr_structure("reference", self.reffile).get_TCRs()
        )[0]
        self.set_reference_residues()

    def set_reffile(self):
        """Super method for setting MHC reference file.

        Raises:
            NotImplementedError
        """
        raise NotImplementedError(
            "TCRCom cannot be insantiated directly, instantiate its subclass"
        )

    def set_reference_residues(self):
        """Set TCR and MHC reference residues."""
        self.set_tcr_reference()
        self.set_mhc_reference()

    def set_tcr_reference(self):
        """Set TCR variable domain residues in reference model as residues numbered 1 to 121 for VA and 1 to 126 for VB."""
        self.reference_VA_residues = [
            r
            for r in self.ref_model.get_VA().get_residues()
            if r.get_id()[1] >= 1 and r.get_id()[1] <= 121
        ]
        self.reference_VB_residues = [
            r
            for r in self.ref_model.get_VB().get_residues()
            if r.get_id()[1] >= 1 and r.get_id()[1] <= 126
        ]

    def set_mhc_reference(self):
        """Super method for setting MHC reference residues for superposition. Overwritten by MHC class sepcific methods.

        Raises:
            NotImplementedError
        """
        raise NotImplementedError(
            "TCRCom cannot be insantiated directly, instantiate its subclass"
        )

    def get_filtered_TCR_residues(self, tcr: "TCR") -> "tuple[list[Bio.PDB.Residue]]":
        """Get variable domain residues of query TCR and filter out those without a counterpart in the reference TCR.

        Args:
            tcr (TCR): TCR structure object

        Returns:
            tuple[list[Bio.PDB.Residue]]: VA_residues, VB_residues.
        """
        tcr_A_residues = [
            tcr.get_VA()[r.get_id()]
            for r in self.reference_VA_residues
            if r.get_id() in tcr.get_VA()
        ]
        tcr_B_residues = [
            tcr.get_VB()[r.get_id()]
            for r in self.reference_VB_residues
            if r.get_id() in tcr.get_VB()
        ]
        return tcr_A_residues, tcr_B_residues

    def center_of_mass(
        self,
        entity: "Union[Bio.PDB.Entity.Entity, list[Bio.PDB.Atom.Atom]]",
        geometric: bool = False,
    ) -> np.array:
        """Calculate the mass weighted or purely geometric centre of mass of an entity or a list of atoms.

        Args:
            entity (Union[Bio.PDB.Entity.Entity, list[Bio.PDB.Atom.Atom]]): Structure entity whose centre of mass will be calculated
            geometric (bool, optional): Whether to calculate geometric mean. Defaults to False.

        Raises:
            ValueError: Checks input type to ensure atoms can be retrieved.
            ValueError: Unknown atoms.

        Returns:
            np.array: centre of mass
        """
        # Structure, Model, Chain, Residue
        if isinstance(entity, Bio.PDB.Entity.Entity):
            atom_list = entity.get_atoms()
        # List of Atoms
        elif hasattr(entity, "__iter__") and [x for x in entity if x.level == "A"]:
            atom_list = entity
        else:
            raise ValueError(
                "Center of Mass can only be calculated from the following objects:\n"
                "Structure, Model, Chain, Residue, list of Atoms."
            )

        masses, positions = zip(*[(atom.mass, atom.coord) for atom in atom_list])
        positions = np.asarray(positions)
        if "ukn" in set(masses) and not geometric:
            raise ValueError(
                "Some Atoms don't have an element assigned.\n"
                "Try adding them manually or calculate the geometrical center of mass instead"
            )

        if geometric:
            return positions.sum(axis=0) / len(atom_list)
        else:
            return np.matmul(np.asarray(masses), positions) / len(atom_list)

    def add_com(
        self,
        mhc_com: np.array,
        tcr_com: np.array,
        VA_com: np.array,
        VB_com: np.array,
        tcr: "TCR",
    ) -> "TCR":
        """
        Function to add pseudoatoms at MHC-CoM, TCR-CoM, and XYZ axis to the output PDB file

        Args:
            mhc_com (np.array): MHC centre of mass
            tcr_com (np.array): TCR centre of mass
            VA_com (np.array): Alpha chain centre of mass
            VB_com (np.array): Beta chain entre of mass
            tcr (TCR): TCR structure object

        Returns:
            TCR: Copy of the original TCR strucutre object with added pseudo-atoms.
        """
        new_structure = tcr.copy()

        # mhc com
        mhc_com_chain = "X"
        new_structure.add(Bio.PDB.Chain.Chain(mhc_com_chain))
        res_id = (" ", 1, " ")
        new_residue = Bio.PDB.Residue.Residue(res_id, "MCM", " ")
        new_atom = Bio.PDB.Atom.Atom("C", mhc_com, 0, 0.0, " ", "C", 1, "C")
        new_residue.add(new_atom)
        new_structure.child_dict[mhc_com_chain].add(new_residue)

        # tcr com
        tcr_com_chain = "Y"
        new_structure.add(Bio.PDB.Chain.Chain(tcr_com_chain))
        pseudo_atom_ids = ["TCM", "ACM", "BCM"]
        tcr_com_list = [tcr_com, VA_com, VB_com]
        for i in range(3):
            res_id = (" ", i + 1, " ")
            new_residue = Bio.PDB.Residue.Residue(res_id, pseudo_atom_ids[i], " ")
            new_atom = Bio.PDB.Atom.Atom("C", tcr_com_list[i], 0, 0.0, " ", "C", 1, "C")
            new_residue.add(new_atom)
            new_structure.child_dict[tcr_com_chain].add(new_residue)

        # X,Y,Z atoms
        pos = [[50, 0, 0], [0, 50, 0], [0, 0, 50]]
        resn = ["X", "Y", "Z"]
        xyz_chain = "Z"
        new_structure.add(Bio.PDB.Chain.Chain(xyz_chain))
        for i in [0, 1, 2]:
            res_id = (" ", i + 1, " ")
            new_residue = Bio.PDB.Residue.Residue(res_id, resn[i], " ")
            new_atom = Bio.PDB.Atom.Atom("O", pos[i], 0, 0.0, " ", "O", 1, "O")
            new_residue.add(new_atom)
            new_structure.child_dict[xyz_chain].add(new_residue)

        return new_structure

    def calculate_centres_of_mass(
        self,
        tcr: abTCR,
        save_aligned_as: str = None,
    ) -> tuple[np.array]:
        """Calculate the TCR and MHC centres of mass of an stcrpy TCR structure object.

        Args:
            tcr (abTCR): TCR structure object
            save_aligned_as (str): Path to same alignment to. If None or False alignment is not saved. Defaults to None.

        Raises:
            NotImplementedError: Alpha Beta TCR compatible only, Gamma Delta TCRs not implemented.

        Returns:
            tuple[np.array]: tcr_com, mhc_com, tcr_VA_com, tcr_VB_com
        """

        assert len(tcr.get_MHC()) > 0, "No MHC associated with TCR"
        if not isinstance(tcr, abTCR):
            raise NotImplementedError(
                f"TCR MHC geometry only implemented for abTCR types, not {type(tcr)}"
            )

        mhc_atoms = [res["CA"] for res in self.get_filtered_MHC_residues(tcr)]
        ref_mhc_atoms = [
            res["CA"]
            for res in self.reference_MHC_residues
            if (res.parent.chain_type, res.get_id())
            in [(a.parent.parent.chain_type, a.parent.get_id()) for a in mhc_atoms]
        ]

        superimposer = Bio.PDB.Superimposer()
        superimposer.set_atoms(ref_mhc_atoms, mhc_atoms)
        self.mhc_alignment_transform = (x.astype("f") for x in superimposer.rotran)
        superimposer.apply(tcr.parent.get_atoms())

        mhc_com = self.center_of_mass(mhc_atoms, geometric=True)

        tcr_VA_residues, tcr_VB_residues = self.get_filtered_TCR_residues(tcr)
        tcr_VA_atoms = [res["CA"] for res in tcr_VA_residues]
        tcr_VA_com = self.center_of_mass(tcr_VA_atoms, geometric=True)

        tcr_VB_atoms = [res["CA"] for res in tcr_VB_residues]
        tcr_VB_com = self.center_of_mass(tcr_VB_atoms, geometric=True)

        tcr_com = self.center_of_mass(tcr_VA_atoms + tcr_VB_atoms, geometric=True)

        if save_aligned_as:
            if not isinstance(save_aligned_as, str):
                save_aligned_as = os.path.join(
                    os.getcwd(), f"{tcr.parent.parent.id}_{tcr.id}.pdb"
                )
            aligned_tcr = self.add_com(mhc_com, tcr_com, tcr_VA_com, tcr_VB_com, tcr)
            io = TCRIO()
            io.save(aligned_tcr, save_as=save_aligned_as)

        # com_distance = tcr_com - mhc_com
        # r = np.sqrt(np.sum(np.square(com_distance)))
        # theta = np.degrees(np.arctan2(com_distance[1], com_distance[0]))
        # phi = np.degrees(np.arccos(com_distance[2] / r))

        # return r, theta, phi

        return tcr_com, mhc_com, tcr_VA_com, tcr_VB_com


class MHCI_TCRCoM(TCRCoM):
    def __init__(self):
        """TCRCoM module for MHC Class I complexes."""
        super().__init__()

    def set_reffile(self):
        """Sets reference file for MHC class I structures"""
        self.reffile = os.path.join(
            os.path.dirname(__file__),
            "reference_data/dock_reference_1_imgt_numbered.pdb",
        )

    def set_mhc_reference(self):
        """Set class I reference MHC residues"""
        mhc = self.ref_model.get_MHC()[0].get_MH1()
        self.reference_MHC_residues = [
            r for r in mhc.get_residues() if r.get_id()[1] >= 1 and r.get_id()[1] <= 179
        ]

    def get_filtered_MHC_residues(self, tcr: "TCR") -> list[Bio.PDB.Residue]:
        """Retrieve MHC residues from query TCR and filter out those whose counterpart is not found in reference.

        Args:
            tcr (TCR): TCR structure object associated with MHC Class I.

        Returns:
            list[Bio.PDB.Residue]: filtered_MHC_residues
        """
        mhc = tcr.get_MHC()[0]
        if not isinstance(mhc, MHCchain):  # handle single MHC chain case
            mhc = mhc.get_MH1()
        filtered_MHC_residues = [
            mhc[ref_res.get_id()]
            for ref_res in self.reference_MHC_residues
            if ref_res.get_id() in mhc
        ]
        return filtered_MHC_residues


class MHCII_TCRCoM(TCRCoM):
    def __init__(self):
        """TCRCoM module for MHC Class II complexes."""
        super().__init__()

    def set_reffile(self):
        """Sets reference file for MHC class II structures"""
        self.reffile = os.path.join(
            os.path.dirname(__file__),
            "reference_data/dock_reference_2_imgt_numbered.pdb",
        )

    def set_mhc_reference(self):
        """Set class II reference MHC residues"""
        mhc = self.ref_model.get_MHC()[0]
        self.reference_MHC_residues = [
            r
            for r in mhc.get_GA().get_residues()
            if r.get_id()[1] >= 1 and r.get_id()[1] <= 88
        ]
        self.reference_MHC_residues.extend(
            [
                r
                for r in mhc.get_GB().get_residues()
                if r.get_id()[1] >= 1 and r.get_id()[1] <= 87
            ]
        )

    def get_filtered_MHC_residues(self, tcr: "TCR") -> list[Bio.PDB.Residue]:
        """Retrieve MHC residues from query TCR and filter out those whose counterpart is not found in reference.

        Args:
            tcr (TCR): TCR structure object associated with MHC Class II.

        Returns:
            list[Bio.PDB.Residue]: filtered_MHC_residues
        """
        mhc = tcr.get_MHC()[0]
        if hasattr(mhc, "get_GA"):
            filtered_MHC_residues = [
                mhc.get_GA()[ref_res.get_id()]
                for ref_res in self.reference_MHC_residues
                if ref_res.parent.chain_type == "GA"
                and ref_res.get_id() in mhc.get_GA()
            ]
        else:
            warnings.warn(f"No GA chain found for MHC class II: {mhc}")
        if hasattr(mhc, "get_GA"):
            filtered_MHC_residues.extend(
                [
                    mhc.get_GB()[ref_res.get_id()]
                    for ref_res in self.reference_MHC_residues
                    if ref_res.parent.chain_type == "GB"
                    and ref_res.get_id() in mhc.get_GB()
                ]
            )
        else:
            warnings.warn(f"No GB chain found for MHC class II: {mhc}")
        return filtered_MHC_residues
