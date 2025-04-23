import os
import re
import warnings
import numpy as np

import Bio
from Bio.PDB.Superimposer import Superimposer

from .. import tcr_processing


class HADDOCKFormatter:

    def __init__(self, save_dir: str = None):
        """Constructor HADDOCK formatting object.

        Args:
            save_dir (str, optional): Path to save formatted files to. Defaults to None.
        """
        self.save_dir = save_dir if save_dir is not None else "."

    def tcr_to_haddock(self, tcr: "TCR"):
        """Bound reformatting of TCR structure object to HADDOCK compatible PDB file.

        Args:
            tcr (TCR): TCR structure object
        """
        self.write_TCR_pdb_file(tcr, self.save_dir)

    def pMHC_to_haddock(self, mhc: "MHC", antigen: list["Antigen"]):
        """Bound reformatting of MHC and antigen structures object to HADDOCK compatible PDB file.

        Args:
            mhc (MHC): MHC structure object
            antigen (Antigen): Antigen structure object
        """
        self.write_antigen_pdb_file(mhc, antigen, self.save_dir)

    def write_TCR_pdb_file(self, tcr: "TCR", save_dir: str):
        """
        Writes TCR structure to a PDB file in a format HADDOCK can deal with.
        Generates a PDB file, a mapping from the old to the new numbering,
            and a list of active residues to restrain the HADDOCK simulation.

        Args:
            tcr (TCR): The TCR structure.
            save_dir (str): The directory to save the files (default is current directory).
        """
        tcr_id = f"{tcr.parent.parent.id}_{tcr.id}"
        new_tcr_structure = Bio.PDB.Model.Model(id=0)
        residue_conversion = {}
        for i, chain in enumerate(tcr.get_chains()):
            residue_conversion[chain.id] = {}
            new_chain = Bio.PDB.Chain.Chain(id=chain.id)
            selected_residues = [
                res
                for res in Bio.PDB.Selection.unfold_entities(chain, "R")
                if res.id[1] in list(range(1, 130))
            ]
            for residue in selected_residues:
                # handle insertion numbering for HADDOCK
                new_residue = residue.copy()
                if new_residue.id[-1] != " ":
                    new_residue.id = (
                        new_residue.id[0],
                        10 * new_residue.id[1]
                        + (200 * i)
                        + imgt_insertion_char_to_int(new_residue.id[-1]),
                        " ",
                    )
                else:
                    new_residue.id = (
                        new_residue.id[0],
                        new_residue.id[1] + (200 * i),
                        new_residue.id[-1],
                    )
                if new_residue.id != residue.id:
                    residue_conversion[chain.id][residue.id] = new_residue.id
                new_chain.add(new_residue)
            new_tcr_structure.add(new_chain)
        if not os.path.exists(os.path.join(save_dir, tcr_id)):
            os.mkdir(os.path.join(save_dir, tcr_id))
        with open(
            os.path.join(save_dir, f"{tcr_id}/{tcr_id}_haddock_active_residues.txt"),
            "w",
        ) as f:
            # get cdr numbering
            active_residues = []
            for chain in tcr.get_chains():
                cdrs = chain.get_CDRs()

                res_list = [r.id for cdr in cdrs for r in cdr.get_residues()]
                for res_key in res_list:
                    # res_key = (" ", *res[0])
                    if res_key in residue_conversion[chain.id]:
                        active_residues.append(residue_conversion[chain.id][res_key][1])
                    else:
                        active_residues.append(res_key[1])
            f.write("TCR ACTIVE RESIDUES FOR HADDOCK\n")
            f.write(",".join([str(r) for r in active_residues]))
            f.write("\n")

        with open(
            os.path.join(save_dir, f"{tcr_id}/{tcr_id}_haddock_renumbering.txt"), "w"
        ) as f:
            f.write("TCR RESIDUE RENUMBERING FOR HADDOCK\n")
            for chain_id in residue_conversion:
                for res, new_res in residue_conversion[chain_id].items():
                    r_as_str = f"({chain_id},({res[0]},{res[1]},{res[2]}),({new_res[0]},{new_res[1]},{new_res[2]})\n"
                    f.write(r_as_str)

        pdb_io = Bio.PDB.PDBIO()
        pdb_io.set_structure(new_tcr_structure)
        filename = os.path.join(save_dir, f"{tcr_id}/{tcr_id}_tcr_for_docking.pdb")
        pdb_io.save(filename)
        return filename

    def write_antigen_pdb_file(
        self, mhc: "MHC", antigen: list["Antigen"], save_dir: str
    ):
        """
        Writes the antigen PDB file for docking with HADDOCK.
        Generates a PDB file, a file containing the renumbering mapping, and a list of active residues to restrict the simulation.

        Args:
            mhc (MHC): MHC structure object.
            antigen (list[Antigen]): List containing antigen chain. Should be length 1.
            save_dir (str, optional): The directory to save the PDB file. Defaults to ".".

        Returns:
            str: The filename of the saved antigen PDB file.
        """
        # structure = p.get_structure()
        # chains = [i for i in structure.get_antigens() if i.id in antigen_chain_ids]
        mhc_chains = [c for c in mhc.get_chains() if c.chain_type != "B2M"]
        antigen_chains = mhc_chains + antigen
        mhc_id = f"{mhc.parent.parent.id}_MHC_{''.join([c.id for c in antigen_chains])}"

        new_antigen_structure = Bio.PDB.Model.Model(id=0)
        residue_conversion = {}
        for i, chain in enumerate(antigen_chains):
            residue_conversion[chain.id] = {}
            new_chain = Bio.PDB.Chain.Chain(id=chain.id)
            for residue in chain.get_residues():
                # handle insertion numbering for HADDOCK
                new_residue = residue.copy()
                new_residue.id = (
                    new_residue.id[0],
                    new_residue.id[1] + (500 * i),
                    new_residue.id[-1],
                )
                if new_residue.id != residue.id:
                    residue_conversion[chain.id][residue.id] = new_residue.id
                new_chain.add(new_residue)
            new_antigen_structure.add(new_chain)

        if not os.path.exists(os.path.join(save_dir, mhc_id)):
            os.mkdir(os.path.join(save_dir, mhc_id))

        with open(
            os.path.join(save_dir, f"{mhc_id}/{mhc_id}_haddock_active_residues.txt"),
            "a",
        ) as f:
            # get peptide numbering and select as active
            active_residues = []
            for chain in antigen:
                res_list = list(r.id for r in chain.get_residues())
                for res_key in res_list:
                    if res_key in residue_conversion[chain.id]:
                        active_residues.append(residue_conversion[chain.id][res_key][1])
                    else:
                        active_residues.append(res_key[1])
            f.write("ANTIGEN ACTIVE RESIDUES FOR HADDOCK\n")
            f.write(",".join([str(r) for r in active_residues]))
            f.write("\n")

        with open(
            os.path.join(save_dir, f"{mhc_id}/{mhc_id}_haddock_renumbering.txt"), "a"
        ) as f:
            f.write("ANTIGEN RESIDUE RENUMBERING FOR HADDOCK\n")
            for chain in residue_conversion:
                for res, new_res in residue_conversion[chain].items():
                    r_as_str = f"({chain},({res[0]},{res[1]},{res[2]}),({new_res[0]},{new_res[1]},{new_res[2]})\n"
                    f.write(r_as_str)

        pdb_io = Bio.PDB.PDBIO()
        pdb_io.set_structure(new_antigen_structure)
        filename = os.path.join(save_dir, f"{mhc_id}/{mhc_id}_antigen_for_docking.pdb")
        pdb_io.save(filename)
        return filename


class HADDOCKResultsParser:

    def __init__(
        self,
        haddock_results_dir: str,
        tcr_renumbering_file: str = None,
        pmhc_renumbering_file: str = None,
    ):
        """Parser for results from HADDOCK simulations. Renumbers TCR, MHC and Antigen using renumbering files, and parses result metrics.

        Args:
            haddock_results_dir (str): path to HADDOCK simulation results.
            tcr_renumbering_file (str, optional): path to text file containing TCR renumbering to restore from HADDOCK compatible numbering. Defaults to None.
            pmhc_renumbering_file (str, optional): path to text file containing MHC and antigen renumbering to restore from HADDOCK compatible numbering. Defaults to None.
        """

        self.haddock_results_dir = haddock_results_dir
        self.tcr_renumbering_file = tcr_renumbering_file
        self.pmhc_renumbering_file = pmhc_renumbering_file

        if self.haddock_results_dir.endswith(".tgz"):
            warnings.warn(
                "HADDOCK results are compressed. Decompress results before proceeding."
            )

    def renumber_all_haddock_predictions(self):
        """Renumber all haddock predictions contained in results folder. Requires standard HADDOCK output directory format."""
        path = os.path.join(self.haddock_results_dir, "structures/it1/")
        pattern = re.compile(r"complex_.*\.pdb")

        for filename in os.listdir(path):
            if pattern.match(filename):
                file_path = os.path.join(path, filename)
                self.renumber_haddock_prediction(
                    file_path,
                    self.tcr_renumbering_file,
                    self.pmhc_renumbering_file,
                )

    def renumber_haddock_prediction(
        self,
        docked_prediction_file: str,
        haddock_renumbering_file: str,
        antigen_renumbering_file: str = None,
    ) -> Bio.PDB.Model.Model:
        """
        Renumber the HADDOCK prediction based on the renumbering files.

        Args:
            docked_prediction_file (str): Path to the docked prediction file.
            haddock_renumbering_file (str): Path to the HADDOCK renumbering file.
            antigen_renumbering_file (str, optional): Path to the antigen renumbering file.
                                                    Needed for TCR only PDBs with no antigen. Defaults to None.

        Returns:
            Bio.PDB.Model.Model: The renumbered HADDOCK prediction.

        Raises:
            ValueError: If the renumbering index is not found in the renumbering file.

        """

        # initialise file parsers
        tcr_parser = tcr_processing.TCRParser.TCRParser()
        bio_parser = Bio.PDB.PDBParser()

        # find chain ID of TCR to distinguish TCR from antigen
        tcr_chain_id = list(
            tcr_parser.get_tcr_structure("tmp", docked_prediction_file).get_TCRchains()
        )[0].get_id()
        docked_prediction = bio_parser.get_structure("docked", docked_prediction_file)

        # get chains of HADDOCK dock
        merged_tcr_chain = docked_prediction[0][tcr_chain_id]
        merged_antigen_chain = [
            chain
            for chain in docked_prediction.get_chains()
            if chain.id != merged_tcr_chain.id
        ][0]

        # get renumbering
        with open(haddock_renumbering_file, "r") as f:
            lines = f.readlines()

        try:
            antigen_renumbering_index = lines.index(
                "ANTIGEN RESIDUE RENUMBERING FOR HADDOCK\n"
            )
            antigen_renumber_indices = (
                antigen_renumbering_index + 1,
                -1,
            )
        except ValueError:
            antigen_renumbering_index = -1
            antigen_renumber_indices = (
                -1,
                -1,
            )
        tcr_renumber_indices = (1, antigen_renumbering_index)

        # if antigen renumbering file is provided, get antigen renumbering from there
        if antigen_renumbering_file is not None:
            lines = (
                lines[: antigen_renumber_indices[0] - 1]
                if antigen_renumber_indices[0] != -1
                else lines
            )
            tcr_renumber_indices = (1, len(lines) - 1)
            antigen_renumber_indices = (len(lines) + 1, -1)
            with open(antigen_renumbering_file, "r") as f:
                antigen_xtal_lines = f.readlines()
            antigen_renumbering_index = antigen_xtal_lines.index(
                "ANTIGEN RESIDUE RENUMBERING FOR HADDOCK\n"
            )
            lines.extend(antigen_xtal_lines[antigen_renumbering_index:])

        # renumber TCR by creating new PDB model and populating with residues
        tcr_parsed_lines = list(
            map(
                parse_renumbered_line,
                lines[tcr_renumber_indices[0] : tcr_renumber_indices[1]],
            )
        )
        changed_tcr_chain_ids, _, _ = list(zip(*tcr_parsed_lines))

        tcr = Bio.PDB.Model.Model(id=0)

        if len(set(changed_tcr_chain_ids)) > 1:
            id_for_conserved_numbered_chain = min(
                set(changed_tcr_chain_ids), key=changed_tcr_chain_ids.count
            )
        else:
            id_for_conserved_numbered_chain = merged_tcr_chain.id
        tcr.add(Bio.PDB.Chain.Chain(id=id_for_conserved_numbered_chain))
        try:
            tcr.add(
                Bio.PDB.Chain.Chain(
                    id=max(set(changed_tcr_chain_ids), key=changed_tcr_chain_ids.count)
                )
            )
            second_tcr_chain_id = None
        except Bio.PDB.PDBExceptions.PDBConstructionException:
            for id_to_try in "ABCDEFGH":
                try:
                    tcr.add(Bio.PDB.Chain.Chain(id=id_to_try))
                    second_tcr_chain_id = id_to_try
                    break
                except Bio.PDB.PDBExceptions.PDBConstructionException:
                    continue

        renumbered_residues = {}
        for renumbering in tcr_parsed_lines:
            try:
                residue = merged_tcr_chain[renumbering[-1]]
                merged_tcr_chain.detach_child(renumbering[-1])
                residue.id = renumbering[1]
                if second_tcr_chain_id is None:
                    if renumbering[0] not in renumbered_residues:
                        renumbered_residues[renumbering[0]] = []
                    renumbered_residues[renumbering[0]].append(residue)
                else:
                    if second_tcr_chain_id not in renumbered_residues:
                        renumbered_residues[second_tcr_chain_id] = []
                    renumbered_residues[second_tcr_chain_id].append(residue)
                # tcr[renumbering[0]].add(residue)
            except KeyError as e:
                warnings.warn(
                    f"""Renumbering {renumbering} failed with Key Error {e}"""
                )
        for residue in merged_tcr_chain.get_residues():
            if id_for_conserved_numbered_chain not in renumbered_residues:
                renumbered_residues[id_for_conserved_numbered_chain] = []
            renumbered_residues[id_for_conserved_numbered_chain].append(residue)
            # tcr[id_for_conserved_numbered_chain].add(residue)

        # sort the residues
        for chain_id in renumbered_residues:
            sorted_residues = sort_residues_by_imgt_numbering(
                renumbered_residues[chain_id]
            )
            for res in sorted_residues:
                tcr[chain_id].add(res)

        # renumber antigen
        antigen_parsed_lines = list(
            map(
                parse_renumbered_line,
                lines[antigen_renumber_indices[0] :],
            )
        )
        changed_antigen_chain_ids, _, _ = list(zip(*antigen_parsed_lines))
        try:
            tcr.add(Bio.PDB.Chain.Chain(id=merged_antigen_chain.id))
        except Bio.PDB.PDBExceptions.PDBConstructionException:
            for id_to_try in "ABCDEFGH":
                if id_to_try == set(changed_antigen_chain_ids).pop():
                    continue
                try:
                    tcr.add(Bio.PDB.Chain.Chain(id=id_to_try))
                    merged_antigen_chain.id = id_to_try
                    break
                except Bio.PDB.PDBExceptions.PDBConstructionException:
                    continue
        assert (
            len(set(changed_antigen_chain_ids)) == 1
        ), "More than one chain renumbered in renumbering file"
        try:
            tcr.add(Bio.PDB.Chain.Chain(id=set(changed_antigen_chain_ids).pop()))
            new_antigen_chain_id = None
        except Bio.PDB.PDBExceptions.PDBConstructionException:
            for id_to_try in "ABCDEFGH":
                if id_to_try == set(changed_antigen_chain_ids).pop():
                    continue
                try:
                    tcr.add(Bio.PDB.Chain.Chain(id=id_to_try))
                    new_antigen_chain_id = id_to_try
                    break
                except Bio.PDB.PDBExceptions.PDBConstructionException:
                    continue

        for renumbering in antigen_parsed_lines:
            try:
                residue = merged_antigen_chain[renumbering[-1]]
                merged_antigen_chain.detach_child(renumbering[-1])
                residue.id = renumbering[1]
                if new_antigen_chain_id is None:
                    tcr[renumbering[0]].add(residue)
                else:
                    tcr[new_antigen_chain_id].add(residue)
            except KeyError as e:
                warnings.warn(
                    f"""Renumbering {renumbering} failed with Key Error {e}"""
                )
        for residue in merged_antigen_chain.get_residues():
            tcr[merged_antigen_chain.id].add(residue)

        # create structure object and save
        tcr_struct = Bio.PDB.Structure.Structure(id=0)
        tcr_struct.add(tcr)

        pdb_io = Bio.PDB.PDBIO()
        pdb_io.set_structure(tcr_struct)
        save_to = "renumbered_" + docked_prediction_file.split("/")[-1]
        filename = os.path.join(*docked_prediction_file.split("/")[:-1], save_to)
        pdb_io.save(filename)

    def get_haddock_scores(self) -> "pandas.DataFrame":
        """Retrieve HADDOCK energy scoes and RMSD evaluations from simulation output:
            \nColumns:
            \n    "haddock_score",
            \n    "interface_rmsd",
            \n    "ligand_rmsd",
            \n    "frac_common_contacts",
            \n    "E_vdw",
            \n    "E_elec",
            \n    "E_air",
            \n    "E_desolv",
            \n    "ligand_rmsd_2",
            \n    "cluster_id",
        Raises:
            FileNotFoundError: HADDOCK file contianing scores not found.

        Returns:
            pandas.DataFrame: DataFrame with HADDOCK simulation metrics.
        """
        import pandas as pd
        import os

        haddock_columns = [
            # 'idx',
            "haddock_score",
            "interface_rmsd",
            "ligand_rmsd",
            "frac_common_contacts",
            "E_vdw",
            "E_elec",
            "E_air",
            "E_desolv",
            "ligand_rmsd_2",
            "cluster_id",
        ]
        haddock_scores_file = "complex_HS_irmsd_lrmsd_fnat.list"
        try:
            df = pd.read_csv(
                os.path.join(self.haddock_results_dir, haddock_scores_file),
                sep=" ",
                names=haddock_columns,
            )
            return df

        except FileNotFoundError:
            raise FileNotFoundError(
                f"File: complex_HS_irmsd_lrmsd_fnat.list containing HADDOCK docking metrics not found in {self.haddock_results_dir}"
            )


def imgt_insertion_char_to_int(char: str) -> int:
    """
    Converts an IMGT insertion character to an integer.

    Args:
        char (str): The IMGT insertion character.

    Returns:
        int: The corresponding integer value.
    """
    return ord(char) - ord("A") + 1


def parse_renumbered_line(line: str) -> tuple:
    """
    Parses a renumbered line from a file and extracts the chain ID, original numbering, and HADDOCK numbering.

    Args:
        line (str): The renumbered line to parse.

    Returns:
        tuple: A tuple containing the chain ID, original numbering, and HADDOCK numbering.

    Example:
        line = "(O,( ,3, ),( ,203, )"
        result = parse_renumbered_line(line)
        # Output: (O)', ('', '3', ''), ('', '203', ''))
    """
    chain_id = line[1]
    content = re.findall(r"\((.*?)\)", line)
    original_numbering = tuple(
        int(x.strip()) if x.isdigit() else x.strip()
        for x in content[0].split("(")[-1].split(",")
    )
    haddock_numbering = tuple(
        int(x.strip()) if x.isdigit() else x.strip()
        for x in re.split(r",\s*", content[1])
    )

    def add_empty_id(numbering):
        return tuple(x if x != "" else " " for x in numbering)

    return chain_id, add_empty_id(original_numbering), add_empty_id(haddock_numbering)


def sort_residues_by_imgt_numbering(
    residues: "list[Bio.PDB.Residue]",
) -> "list[Bio.PDB.Residue]":
    """Sort residues in order by IMGT numbering.

    Args:
        residues (list[Bio.PDB.Residue]): List of IMGT numbered residues.

    Returns:
        list[Bio.PDB.Residue]: Sorted list of IMGT numbered residuess.
    """
    sorted_residues = sorted(residues, key=lambda x: (x.id[1], x.id[2]))
    imgt_nr_112_subsequence = [
        (i, res) for i, res in enumerate(sorted_residues) if res.id[1] == 112
    ]
    if len(imgt_nr_112_subsequence) > 0:
        indices, imgt_nr_112_subsequence = list(zip(*imgt_nr_112_subsequence))
        sorted_imgt_nr_112_subsequence = sorted(
            imgt_nr_112_subsequence, key=lambda x: x.id[2], reverse=True
        )
        for i, idx in enumerate(indices):
            sorted_residues[idx] = sorted_imgt_nr_112_subsequence[i]
    return sorted_residues
