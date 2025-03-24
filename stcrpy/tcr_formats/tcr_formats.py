from Bio.SeqUtils import seq1
import json
import os


def to_AF3_json(
    tcr: "TCR",
    tcr_only: bool = True,
    save: bool = True,
    save_dir: str = "",
    name: str = None,
    V_domain_only: bool = False,
) -> dict:
    """Converts TCR object to dict in Alphafold 3 JSON input format, ie. amino acid sequences.
    Eg:
    {
        "name": Job name,
        "modelSeeds": [],
        "sequences": [
            {"proteinChain": {"sequence": AAAAAAAAAAAAAA, "count": 1}},
            {"proteinChain": {"sequence": AAAAAAAAAAAAAA, "count": 1}},
            {"proteinChain": {"sequence": AAAAAAAAAAAAAA, "count": 1}},
        ],
    }

    Args:
        tcr (TCR): TCR structure object
        tcr_only (bool, optional): Whether to include TCR sequence only, excluding antigen and MHC. Defaults to True.
        save (bool, optional): Whether to save dict as JSON file. Defaults to True.
        save_dir (str, optional): Directory to save JSON files to. Defaults to "".
        name (str, optional): TCR ID to use as name for AF3 job. Defaults to None.
        V_domain_only (bool, optional): Include full TCR sequence or only the variable domain (1-128 IMGT numbering). Defaults to False.

    Returns:
        dict: Nested dictionary of AF3 sequence inputs.
    """
    if V_domain_only:
        residue_nrs = list(range(128))
    else:
        residue_nrs = None
    tcr_sequences = get_sequences(tcr, residues_to_include=residue_nrs)
    if not tcr_only:
        if len(tcr.get_MHC()) > 0:
            mhc_sequences = get_sequences(tcr.get_MHC()[0])
            tcr_sequences.update(mhc_sequences)

        if len(tcr.get_antigen()) > 0:
            antigen_sequence = get_sequences(tcr.get_antigen()[0])
            tcr_sequences.update(antigen_sequence)
    name = name if name is not None else f"{tcr.parent.parent.id}_{tcr.id}"
    tcr_json = {
        "name": name,
        "modelSeeds": [],
        "sequences": [
            {"proteinChain": {"sequence": seq, "count": 1}}
            for _, seq in tcr_sequences.items()
        ],
    }
    if save:
        path = os.path.join(save_dir, f"{name}.json")
        with open(path, "w") as f:
            json.dump(tcr_json, f)
    return tcr_json


def get_sequences(
    entity: "Bio.PDB.Entity",
    amino_acids_only: bool = True,
    residues_to_include: list = None,
) -> dict:
    """Extract seqeunces from strcuture objects as dictionary.

    Args:
        entity (Bio.PDB.Entity): Stucture object
        amino_acids_only (bool, optional): Whether to remove non-amino acid 'X' from sequences. Defaults to True.
        residues_to_include (list, optional): List of residue IDs to include in sequence. Defaults to None.

    Raises:
        e: AttributeError if entity has no attribute .get_chains(). The assuems entity is chain level and returns single sequence

    Returns:
        dict: Dictionary of amino acid sequences, keyed by chain ID in strcuctre entity.
    """

    if residues_to_include is None:

        def residue_filter(res):
            return True

    else:

        def residue_filter(res):
            return res.id[1] in residues_to_include
    try:
        sequences = {
            chain.id: seq1(
                "".join(residue.resname for residue in chain if residue_filter(residue))
            )
            for chain in entity.get_chains()
        }
    except AttributeError as e:
        if entity.level == "C":
            sequences = {
                entity.id: seq1(
                    "".join(
                        residue.resname for residue in entity if residue_filter(residue)
                    )
                )
            }
        else:
            raise e
    if amino_acids_only:
        sequences = {k: seq.replace("X", "") for k, seq in sequences.items()}
    return sequences
