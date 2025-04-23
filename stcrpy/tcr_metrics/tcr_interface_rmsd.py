import warnings
import Bio
from Bio.PDB.Superimposer import Superimposer
import numpy as np

class InterfaceRMSD:
    def __init__(self):
        return

    def get_interface_rmsd(self, dock: "abTCR", reference: "abTCR") -> float:
        """
        Calculates the root-mean-square deviation (RMSD) between the interface residues of a
        docked TCR structure and a reference TCR structure.

        Args:
            dock (abTCR): The docked TCR structure.
            reference (abTCR): The reference TCR structure.

        Returns:
            float: The RMSD value between the interface residues of the docked TCR structure and
                    the reference TCR structure.
        """

        # check all residues in dock can be mapped to reference
        dock_to_ref_chain_mapping = self.check_residue_mapping(dock, reference)

        # get the interface of the reference
        tcr_interface, antigen_interface = self.get_interface_residues(reference)

        # align the dock by MHC
        self.align_by_mhc(dock, reference, dock_to_ref_chain_mapping)

        # get the docked residues found in the reference interface
        ref_to_dock_chain_mapping = {v: k for k, v in dock_to_ref_chain_mapping.items()}

        try:
            docked_tcr_interface = [
                dock[ref_to_dock_chain_mapping[res.parent.id]][res.id]
                for res in tcr_interface
            ]
            docked_antigen_interface = [
                dock.parent[ref_to_dock_chain_mapping[res.parent.id]][res.id]
                for res in antigen_interface
            ]
        except KeyError as e:
            warnings.warn(
                f"""Key error {str(e)} matching dock chains to reference chains for dock: {
                    ' '.join([str(i) for i in dock.full_id])
                } to reference {
                    ' '.join([str(i) for i in reference.full_id])
                }. Interface RMSD could not be calculated."""
            )
            return None

        # extract coordinates from interfaces
        reference_coordinates = np.asarray(
            [
                atom.get_coord()
                for res in (tcr_interface + antigen_interface)
                for atom in res
                if atom.element in ["N", "O", "C", "S"]
                and (
                    atom.id
                    in [
                        a.id
                        for a in dock.parent[ref_to_dock_chain_mapping[res.parent.id]][
                            res.id
                        ].get_atoms()
                    ]
                )
            ]
        )

        docked_coordinates = np.asarray(
            [
                atom.get_coord()
                for res in (docked_tcr_interface + docked_antigen_interface)
                for atom in res
                if atom.element in ["N", "O", "C", "S"]
                and (
                    atom.id
                    in [
                        a.id
                        for a in reference.parent[
                            dock_to_ref_chain_mapping[res.parent.id]
                        ][res.id].get_atoms()
                    ]
                )
            ]
        )

        # calculate rmsd
        rmsd = np.sqrt(
            ((docked_coordinates - reference_coordinates) ** 2).sum()
            / len(docked_coordinates)
        )
        return rmsd

    def check_residue_mapping(self, dock: "abTCR", reference: "abTCR") -> dict:
        chain_mapping = {}

        for i, tcr_chain in enumerate(dock.get_chains()):
            try:
                for j, res in enumerate(tcr_chain.get_residues()):
                    if j > 3 and j < len(tcr_chain) - 2:
                        assert (
                            res.resname
                            == list(reference.get_chains())[i][res.id].resname
                        ), f"""
    TCR chain mapping {tcr_chain.id} -> {list(reference.get_chains())[i]} failed. Trying chain swap."""
                        chain_mapping[tcr_chain.id] = list(reference.get_chains())[i].id
            except (AssertionError, KeyError):
                for j, res in enumerate(tcr_chain.get_residues()):
                    if (
                        j > 3 and j < len(tcr_chain) - 2
                    ):  # avoids small mismatches at beginnings and ends of sequences
                        assert (
                            res.resname
                            == list(reference.get_chains())[1 - i][res.id].resname
                        ), f"""
    TCR chain mapping {tcr_chain.id} -> {list(reference.get_chains())[1-i]} failed. Residue mapping failed"""
                        chain_mapping[tcr_chain.id] = list(reference.get_chains())[
                            1 - i
                        ].id

        for i, antigen_chain in enumerate(dock.get_antigen()):
            try:
                for res in antigen_chain.get_residues():
                    assert (
                        res.resname == reference.get_antigen()[i][res.id].resname
                    ), f"""
    antigen chain mapping {antigen_chain.id} -> {reference.get_antigen()[i]} failed. Trying chain swap."""
                    chain_mapping[antigen_chain.id] = reference.get_antigen()[i].id
            except (AssertionError, KeyError):
                for res in antigen_chain.get_residues():
                    assert (
                        res.resname == reference.get_antigen()[1 - i][res.id].resname
                    ), f"""
    antigen chain mapping {antigen_chain.id} -> {reference.get_antigen()[1-i]} failed. Residue mapping failed"""
                    chain_mapping[antigen_chain.id] = reference.get_antigen()[1 - i].id
        return chain_mapping

    def get_interface_residues(
        self, tcr: "abTCR", angstrom_cutoff: float = 8.0
    ) -> list:
        """
        Retrieves the interface residues between a TCR and its antigen based on a distance cutoff.

        Args:
            tcr (abTCR): The TCR object.
            angstrom_cutoff (float, optional): The distance cutoff in angstroms. Defaults to 8.0.

        Returns:
            tuple: A tuple containing two lists: the interface residues of the TCR and the
                    interface residues of the antigen.
        """
        tcr_c_alphas = [atom for atom in tcr.get_atoms() if atom.id == "CA"]
        antigen_c_alphas = [
            atom
            for chain in tcr.get_antigen()
            for atom in chain.get_atoms()
            if atom.id == "CA"
        ]

        tcr_c_coords = np.asarray([[x.get_coord()] for x in tcr_c_alphas])
        antigen_c_coords = np.asarray([[x.get_coord() for x in antigen_c_alphas]])

        tcr_c_coords = np.broadcast_to(
            tcr_c_coords, (tcr_c_coords.shape[0], antigen_c_coords.shape[1], 3)
        )

        antigen_c_coords = np.broadcast_to(
            antigen_c_coords, (tcr_c_coords.shape[0], antigen_c_coords.shape[1], 3)
        )

        pairwise_distances = np.sqrt(((tcr_c_coords - antigen_c_coords) ** 2).sum(-1))
        contacts = np.argwhere(pairwise_distances <= angstrom_cutoff)
        tcr_interface_idx = set(contacts[:, 0])
        antigen_interface_idx = set(contacts[:, 1])
        tcr_interface = [tcr_c_alphas[idx].parent for idx in tcr_interface_idx]
        antigen_interface = [
            antigen_c_alphas[idx].parent for idx in antigen_interface_idx
        ]

        return tcr_interface, antigen_interface

    def align_by_mhc(
        self, dock: "abTCR", reference: "abTCR", chain_mapping: dict
    ) -> None:
        """
        Aligns the docked TCR structure to the reference TCR structure by aligning the MHC.

        Args:
            dock (abTCR): The docked TCR structure.
            reference (abTCR): The reference TCR structure.
            chain_mapping (dict): A dictionary mapping the chain IDs of the docked TCR structure to the chain IDs of the
                                    reference TCR structure.

        Returns:
            None
        """
        mhc_chain = dock.get_MHC()
        assert len(mhc_chain) >= 1, ValueError("No MHC chains found")
        if hasattr(mhc_chain[0], "get_MH1"):
            mhc_chain = mhc_chain[
                0
            ].get_MH1()  # This will only work for class I MHC, ie. single chain helices.
            reference_mhc_chain = reference.get_MHC()[0].get_MH1()
        else:
            # For Class II MHC try creating new entity with GA and GB chains
            class_II_mhc_chain = Bio.PDB.Entity.Entity()
            class_II_mhc_chain.add(mhc_chain[0].get_GA())
            class_II_mhc_chain.add(mhc_chain[0].get_GB())
            mhc_chain = class_II_mhc_chain
            reference_mhc_chain = Bio.PDB.Entity.Entity()
            reference_mhc_chain.add(reference.get_MHC()[0].get_GA())
            reference_mhc_chain.add(reference.get_MHC()[0].get_GB())

        mutual_residue_ids = set(
            [r.id for r in reference_mhc_chain.get_residues()]
        ).intersection(set([r.id for r in mhc_chain.get_residues()]))
        reference_atoms = [
            a
            for res in mutual_residue_ids
            for a in reference_mhc_chain[res].get_atoms()
            if a.id in ["N", "C", "O", "CA"]
        ]
        docked_atoms = [
            a
            for res in mutual_residue_ids
            for a in mhc_chain[res].get_atoms()
            if a.id in ["N", "C", "O", "CA"]
        ]

        superimposer = Superimposer()
        superimposer.set_atoms(reference_atoms, docked_atoms)
        superimposer.apply(dock.parent.get_atoms())
