import warnings
import requests
import os

from ..tcr_processing.TCRParser import TCRParser
from .tcr_batch_operations import batch_load_TCRs


def load_TCRs(tcr_structure_files, tcr_ids=None):
    tcr_parser = TCRParser()
    if isinstance(tcr_structure_files, str):  # loading single file
        tcr_id = tcr_structure_files.split("/")[-1].split(".")[
            0
        ]  # set tcr_id to file name without extension
        if tcr_ids is not None:
            if not isinstance(tcr_ids, str):
                warnings.warn(f"TCR ID: {tcr_ids} for a single TCR should be type str.")
            tcr_id = tcr_ids

        tcr_structure = tcr_parser.get_tcr_structure(tcr_id, tcr_structure_files)
        return list(tcr_structure.get_TCRs())

    if len(tcr_structure_files) > 10:
        warnings.warn(
            "Loading more than 10 TCR structure objects into memory. Consider applying generator methods to reduce memory load."
        )

    if tcr_ids is not None:
        if len(tcr_structure_files) == len(tcr_ids):
            return batch_load_TCRs(dict(zip(tcr_ids, tcr_structure_files)))
        else:
            warnings.warn(
                f"Length of TCR IDs {len(tcr_ids)} does not match length of files {len(tcr_structure_files)}. TCR IDs reverted to default."
            )
    return batch_load_TCRs(tcr_structure_files)


def fetch_TCR(pdb_id: str):
    """
    Fetches and parses a T-cell receptor (TCR) structure from the STCRDab or RCSB PDB databases.

    The function first attempts to download a PDB file from the STCRDab database.
    If the PDB file is not found, it falls back to downloading a CIF file from RCSB PDB.
    The downloaded file is then parsed using `TCRParser` to extract TCR structures.

    Parameters:
        pdb_id (str): The PDB identifier of the structure to be fetched.

    Returns:
        - A single TCR structure if exactly one is found.
        - A list of TCR structures if multiple are found.
        - None if no TCRs are identified (with a `UserWarning` issued).

    Raises:
        - A warning if no TCR structures are found in the downloaded file.
        - Prints an error message if the file cannot be downloaded.

    Notes:
        - STCRDab returns an error message if the requested PDB ID does not exist.
        - The function temporarily saves the downloaded file and deletes it after parsing.

    Example:
        tcr = fetch_TCR("6eqa")

    """

    stcrdab_base_url = "https://opig.stats.ox.ac.uk/webapps/stcrdab-stcrpred/pdb/"
    pdb_base_url = "https://files.rcsb.org/download/"

    filename = f"{pdb_id.upper()}.pdb"

    url = stcrdab_base_url + pdb_id.lower()
    response = requests.get(url, stream=True)

    TCR_FOUND = False

    if response.status_code == 200:
        with open(filename, "wb") as file:
            for chunk in response.iter_content(chunk_size=1024):
                file.write(chunk)
            if (
                not b"does not exist" in chunk
            ):  # STCRDab returns '$PDB does not exist for downloading' if PDB code not found in database
                TCR_FOUND = True
    if not TCR_FOUND:
        os.remove(filename)  # remove the file written with response from STCRDab
        filename = f"{pdb_id.upper()}.cif"
        url = pdb_base_url + filename
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(filename, "wb") as file:
                for chunk in response.iter_content(chunk_size=1024):
                    file.write(chunk)
        else:
            print("Failed to download file")

    tcr_parser = TCRParser()
    tcr = list(tcr_parser.get_tcr_structure(pdb_id, filename).get_TCRs())
    os.remove(filename)
    if len(tcr) == 1:
        return tcr[0]
    elif len(tcr) == 0:
        warnings.warn(f"No TCRs identified in {pdb_id}")
        return None
    else:
        return tcr
