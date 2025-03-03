from Bio.SeqUtils import seq1
import json
import os


def to_AF3_json(
    tcr, tcr_only=True, save=True, save_dir="", name=None, V_domain_only=False
):
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


def get_sequences(entity, amino_acids_only=True, residues_to_include=None):

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
