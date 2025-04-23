import warnings
import numpy as np
import pandas as pd

from Bio.PDB.Superimposer import Superimposer

from . import constants


class RMSD:
    def __init__(self):
        return

    @staticmethod
    def _retrieve_chain(tcr, chain_type):
        assert chain_type in ["A", "B", "G", "D"], ValueError(
            "TCR chain type not recognised"
        )
        try:
            return tcr[tcr.get_domain_assignment()[f"V{chain_type}"]]
        except KeyError:
            # map chain type A to G and B to D
            chain_type = {"A": "G", "B": "D"}[chain_type]
            return tcr[tcr.get_domain_assignment()[f"V{chain_type}"]]

    @staticmethod
    def _rmsd(x1, x2):
        assert x1.shape == x2.shape
        assert x1.shape[-1] == 3
        return np.sqrt(np.mean(3 * (x1 - x2) ** 2))

    def calculate_rmsd(self, tcr_to_align, tcr_ref, save_alignment=False):
        rmsds = {}
        for chain_type in ["A", "B"]:
            chain_to_align = self._retrieve_chain(tcr_to_align, chain_type)
            ref_chain = self._retrieve_chain(tcr_ref, chain_type)

            ref_residue_numbering = [
                x.id
                for x in ref_chain.get_residues()
                if all([not a.is_disordered() for a in x.get_atoms()])
            ]
            residue_numbering_intersection = [
                x.id
                for x in chain_to_align.get_residues()
                if x.id in ref_residue_numbering
                and all([not a.is_disordered() for a in x.get_atoms()])
            ]

            # Get residues to align
            ref_residues = [
                ref_chain[x]
                for x in residue_numbering_intersection
                if (x in chain_to_align) and (x in ref_chain)
            ]
            to_align_residues = [
                chain_to_align[x]
                for x in residue_numbering_intersection
                if (x in chain_to_align) and (x in ref_chain)
            ]

            # Get backbone atoms to align
            fixed = []
            moved = []
            for i in range(len(to_align_residues)):
                fixed += [
                    ref_residues[i][atom]
                    for atom in constants.ATOM_TYPES[:4]
                    if (atom in to_align_residues[i]) and (atom in ref_residues[i])
                ]
                moved += [
                    to_align_residues[i][atom]
                    for atom in constants.ATOM_TYPES[:4]
                    if (atom in to_align_residues[i]) and (atom in ref_residues[i])
                ]

            # Calculate superimposer and align
            imposer = Superimposer()
            imposer.set_atoms(fixed, moved)
            imposer.apply(tcr_to_align.get_atoms())

            rmsds[chain_type] = (
                imposer.rms
            )  # whole chain RMSD after alignment calculated across all atoms

            if save_alignment:
                tcr_ref.save(
                    save_as=f"{tcr_ref.parent.parent.id}_RMSD_reference_alignment_{chain_type}.pdb",
                    tcr_only=True,
                )
                tcr_to_align.save(
                    save_as=f"{tcr_to_align.parent.parent.id}_RMSD_aligned_to_{tcr_ref.id}_{chain_type}.pdb",
                    tcr_only=True,
                )

            # calculate CDR RMSD
            for CDR_loop_nr, (ref_CDR, aligned_CDR) in enumerate(
                zip(ref_chain.get_CDRs(), chain_to_align.get_CDRs())
            ):
                ref_CDR_atom_coords = np.asarray(
                    [
                        r[a].get_coord()
                        for r in ref_CDR
                        for a in constants.ATOM_TYPES[:4]
                        if r.id in residue_numbering_intersection
                        and a in aligned_CDR[r.id]
                    ]
                )
                aligned_CDR_atom_coords = np.asarray(
                    [
                        r[a].get_coord()
                        for r in aligned_CDR
                        for a in constants.ATOM_TYPES[:4]
                        if r.id in residue_numbering_intersection and a in ref_CDR[r.id]
                    ]
                )

                rmsds[f"CDR{chain_type}{CDR_loop_nr + 1}"] = self._rmsd(
                    ref_CDR_atom_coords, aligned_CDR_atom_coords
                )

            # calculate framework RMSD
            ref_framework_residues = {
                r.id: r
                for fw in ref_chain.get_frameworks()
                for r in fw.get_residues()
                if r.id in residue_numbering_intersection
            }
            aligned_framework_residues = {
                r.id: r
                for fw in chain_to_align.get_frameworks()
                for r in fw.get_residues()
                if r.id in residue_numbering_intersection
            }

            ref_framework_atom_coords = np.asarray(
                [
                    r[a].get_coord()
                    for r_id, r in ref_framework_residues.items()
                    for a in constants.ATOM_TYPES[:4]
                    if r_id in residue_numbering_intersection
                    and a in aligned_framework_residues[r_id]
                ]
            )
            aligned_framework_atom_coords = np.asarray(
                [
                    r[a].get_coord()
                    for r_id, r in aligned_framework_residues.items()
                    for a in constants.ATOM_TYPES[:4]
                    if r_id in residue_numbering_intersection
                    and a in ref_framework_residues[r_id]
                ]
            )

            rmsds[f"FW{chain_type}"] = self._rmsd(
                ref_framework_atom_coords, aligned_framework_atom_coords
            )
        return rmsds

    def rmsd_from_files(self, pred_and_target_files: list) -> pd.DataFrame:
        """Calculates the RMSD between TCR structures from a list of files.

        Args:
            pred_and_target_files (list of tuples): List of tuples, where each tuple contains
                the path to the predicticted PDB at index 0 and the path to the target PDB at index 1.

        Returns:
            pandas.Dataframe: Pandas dataframe indexed by the TCR ID of the predicted structure, with columns
                containing the RMSD of the whole alpha and beta chain, and all framework and CDR regions.
        """
        from ..tcr_methods.tcr_methods import load_TCRs

        all_rmsds = {}
        for pred_tcr_file, target_tcr_file in pred_and_target_files:
            pred_tcr, target_tcr = load_TCRs([pred_tcr_file, target_tcr_file])
            all_rmsds[pred_tcr.parent.parent.id] = self.calculate_rmsd(
                pred_tcr, target_tcr
            )
        return pd.DataFrame(all_rmsds).transpose()
