"""
Created on 10 May 2017
@author: leem

Implementation to call anarci (built-in to STrDab) to annotate structures.
"""
import sys
import warnings

from Bio.PDB.Polypeptide import aa1, aa3  # to allow me to return "X" if not found.

to_one_letter_code = dict(list(zip(aa3, aa1)))

# Import TCRDB's constants and common functions.
from .utils.constants import TCR_CHAINS

NUMBERING_OFFSET = 1000
'''Offset needed to number things like scTCRs where two numberings are present and they may clash.'''


def call_anarci(
    seq,
    allow=set(
        [
            "B",
            "A",
            "D",
            "G",
            "GA1",
            "GA2",
            "GA1L",
            "GA2L",
            "GA",
            "GB",
            "B2M",
            "MH1",
            "MR1",
            "MR2",
        ]
    ),
) -> tuple[list[tuple[tuple[int, str], str]], list[str], dict[str, list[tuple[str, str], float]]]:
    """
    Use the ANARCI program to number the sequence.

    Args:
        seq: An amino acid sequence that you wish to number.

    Returns:
        numbering, chain type, germline information
    """
    from anarci import number as anarci_number

    finished_numbering = False

    numberings = []
    chain_types = []
    germline_infos = {}

    seq_position = 0
    seq_length = len(seq)

    while not finished_numbering:
        numbering, chain_type, germline_info = anarci_number(
            seq[seq_position:], allow=allow, assign_germline=True
        )

        if numbering and chain_type in allow:
            clean_numbering = [(_, aa) for _, aa in numbering if aa != "-"]

            # Needed for scTCRs
            clean_seq_ids = {seq_id for seq_id, _ in clean_numbering}
            numberings_seq_ids = {seq_id for seq_id, _ in numberings}

            if len(clean_seq_ids & numberings_seq_ids) > 0:
                clean_numbering = [
                    ((res_seq_id + NUMBERING_OFFSET, insert_code), res_olc)
                    for (res_seq_id, insert_code), res_olc in clean_numbering
                ]

            numberings.extend(clean_numbering)
            chain_types.append(chain_type)
            germline_infos.update(germline_info)

            clean_sequence = ''.join([olc for _, olc in clean_numbering])
            seq_position += seq.index(clean_sequence) + len(clean_sequence)

            if seq_position >= seq_length:
                finished_numbering = True

        else:
            finished_numbering = True

    return numberings, chain_types, germline_infos


def annotate(chain):
    """
    Annotate the sequence of a chain object from TCRDB.TcrPDB
    # e.g. if you have chains B, A and X, you want to force the annotator to return the annotation
    # for B and A but not for X (the antigen)

    returns a dictionary which has the residue ids as key and the annotation as value or is False,
    and chain type which is B/A/G/D/MH1/GA/GB/B2M or False.
    """
    sequence_list, sequence_str = extract_sequence(chain)
    numbering, chain_types, germline_info = call_anarci(sequence_str)

    aligned_numbering = align_numbering(numbering, sequence_list)

    if ('A' in chain_types and 'B' in chain_types) or ('D' in chain_types and 'G' in chain_types):
        scTCR = True
        chain_type = ''.join(chain_types)

        domain1 = {
            seq_id_og: (res_seq_id_number, insert_code_number)
            for seq_id_og, (res_seq_id_number, insert_code_number) in aligned_numbering.items()
            if res_seq_id_number < NUMBERING_OFFSET
        }
        domain2 = {
            seq_id_og: (res_seq_id_number - NUMBERING_OFFSET, insert_code_number)
            for seq_id_og, (res_seq_id_number, insert_code_number) in aligned_numbering.items()
            if res_seq_id_number >= NUMBERING_OFFSET
        }

        aligned_numbering = (domain1, domain2)

    else:
        scTCR = False
        chain_type = '/'.join(chain_types)

        if chain_type == 'GA1/GA2':
            chain_type = 'MH1'

        elif chain_type == 'MR1/MR2':
            chain_type = 'MR1'

        elif chain_type == 'GA1L/GA2L':
            chain_type = 'CD1'

    return aligned_numbering, chain_type, germline_info, scTCR


def extract_sequence(
    chain, selection=False, return_warnings=False, ignore_hets=False, backbone=False
):
    """
    Get the amino acid sequence of the chain.
    Residues containing HETATOMs are skipped -->  Residues containing HETATOMs are checked as an amino acid.

    Residues containing HETATOMs are checked  to be amino acids and the single letter returned.

    This works provided the residues in the chain are in the correct order.

    Args:
        selection: a selection object to select certain residues
        return_warnings: Flag to return a list of warnings or not
        backbone: Flag whether to only show residues with a complete backbone (in the structure) or not.
    Returns:
        The sequence in a resid:aa tuple list and the sequence as a string.

    """
    sequence_list = []
    warnings = []
    for residue in chain.get_list():
        if (
            residue.id[0] != " "
        ):  # skip HETATOMs - this is not necesserily a good idea, flag to the user that is has been done.
            #            if residue.get_resname() not in to_one_letter_code: # Check that the residue can be converted into a single letter.
            #                continue
            #            if residue.get_resname() in to_one_letter_code: # Check that the residue can be converted into a single letter.
            #                pass
            if residue.get_resname() in to_one_letter_code:
                if ignore_hets:
                    if return_warnings:
                        warnings.append(
                            """Warning: HETATM residue %s at position %s (PDB numbering) found in chain %s.
                            Not including it in structure's sequence."""
                            % (
                                residue.get_resname(),
                                str(residue.id[1]) + residue.id[2].strip(),
                                residue.parent.id,
                            )
                        )
                    else:
                        sys.stderr.write(
                            """Warning: HETATM residue %s position %s (PDB numbering) found in chain %s.
                            Not including it in structure's sequence.\n"""
                            % (
                                residue.get_resname(),
                                str(residue.id[1]) + residue.id[2].strip(),
                                residue.parent.id,
                            )
                        )
                    continue
            else:
                continue

        if selection:
            if not selection.accept(residue):
                continue

        atoms_of_residue = list(residue.child_dict.keys())
        backboneCondition = (
            "N" in atoms_of_residue
            and "C" in atoms_of_residue
            and "CA" in atoms_of_residue
            and "O" in atoms_of_residue
        )  # Boolean to hold if residue has a full backbone

        # CASE 1: backbone = True, and residue has a full backbone; convert a.a into single letter
        if backbone and backboneCondition:
            sequence_list.append(
                (residue.id, to_one_letter_code.get(residue.get_resname(), "X"))
            )
        # CASE 2: backbone = True, but residue does not have a full backbone; use a gap in sequence annotation
        elif backbone and not backboneCondition:
            sequence_list.append((residue.id, "-"))
        # CASE 0 (default): don't care about backbone, just write it to sequence if it's found in structure.
        elif not backbone:
            sequence_list.append(
                (residue.id, to_one_letter_code.get(residue.get_resname(), "X"))
            )  # i am

    sequence_str = "".join([r[1] for r in sequence_list])
    if not return_warnings:
        return sequence_list, sequence_str
    else:
        return sequence_list, sequence_str, warnings


def interpret(x):
    """
    Function to interpret an annotation in the form H100A into the form ( 100, 'A' )
    """
    assert x[0] in TCR_CHAINS, x
    try:
        return (int(x[1:]), " ")
    except ValueError:
        return (int(x[1:-1]), x[-1])


class AlignmentError(Exception):
    """Raised when there is an error aligning two sequences."""
    def __init__(self, sequence1: str, sequence2: str) -> None:
        super().__init__(f'Could not align sequences: {sequence1} and {sequence2}')


def align_numbering(numbering, sequence_list):
    """
    Align the sequence that has been numbered to the sequence you input.
    The numbered sequence should be "in" the input sequence.
    If not, supply an alignment dictionary.(align sequences and use get_alignment_dict(ali1,ali2))

    Raises:
        AlignmentError: if the two sequences cannot be aligned
    """
    if not numbering:
        return False

    numbered_sequence = "".join([r[1] for r in numbering])
    input_sequence = "".join([r[1] for r in sequence_list])

    try:
        numbered_sequence_ali, input_sequence_ali = pairwise_alignment(
            numbered_sequence, input_sequence
        )
    except Exception:
        raise AlignmentError(numbered_sequence, input_sequence)

    input_index = 0
    numbered_index = 0

    seq_id = (1, ' ')
    start_flag = True
    numbered_positions = {seq_id for seq_id, _ in numbering}

    aligned_annotations = {}
    for numbered_residue, input_residue in zip(numbered_sequence_ali, input_sequence_ali, strict=True):
        if numbered_residue == input_residue:
            seq_id = numbering[numbered_index][0]

            numbered_index += 1
            start_flag = False

        elif numbered_residue == '-':
            propposed_seq_id = seq_id if start_flag else (seq_id[0] + 1, ' ')

            if propposed_seq_id in numbered_positions:
                for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
                    propposed_seq_id = (seq_id[0], letter)

                    if propposed_seq_id not in numbered_positions:
                        seq_id = propposed_seq_id
                        break

                else:
                    raise AlignmentError(numbered_sequence_ali, input_sequence_ali)

            seq_id = propposed_seq_id

        else:
            raise AlignmentError(numbered_sequence_ali, input_sequence_ali)

        aligned_annotations[sequence_list[input_index][0]] = seq_id
        numbered_positions.add(seq_id)
        input_index += 1

    return aligned_annotations


def pairwise_alignment(seq1, seq2, exact=False):
    """
    Function to do alignment of sequences between sequences using biopython.
    """
    with warnings.catch_warnings():  # prevents pairwise2 deprecation warning from being raised
        warnings.simplefilter("ignore")
        from Bio.pairwise2 import align

    alignment = None
    s1_aln, s2_aln = easy_alignment(seq1, seq2)
    if s1_aln:
        return s1_aln, s2_aln

    if exact:
        # Align with a match score of 1, mismatch of 0, gap opening of -1.001, and gap extension of -1
        alignment = align.globalms(seq1, seq2, 1, 0, -1, -1.001)
    else:
        alignment = align.globalxx(seq1, seq2)

    if alignment:
        aligned_seqs = alignment[0]
        return aligned_seqs[0], aligned_seqs[1]
    else:
        return False, False


def easy_alignment(seq1, seq2):
    """
    Function to align two sequences by checking if one is in the other.
    This function will conserve gaps.
    """
    assert (
        type(seq1) is str and type(seq2) is str
    ), "Sequences must be strings for easy_alignment"
    if seq1 in seq2:
        start = seq2.index(seq1)
        seq1_ali = "-" * start + seq1 + "-" * (len(seq2) - start - len(seq1))
        return seq1_ali, seq2

    elif seq2 in seq1:
        start = seq1.index(seq2)
        seq2_ali = "-" * start + seq2 + "-" * (len(seq1) - start - len(seq2))
        return seq1, seq2_ali

    else:
        # Can't align them # I return just one value here.
        return False, False


def validate_sequence(seq):
    """
    Check whether a sequence is a protein sequence or if someone has submitted something nasty.
    """
    if len(seq) > 10000:
        raise AssertionError("Sequence too long.")
    if any([1 for s in seq.upper() if s not in aa1]):
        raise AssertionError(
            "Unknown amino acid letter found in sequence: " + seq.upper()
        )
    else:
        return True


if __name__ == "__main__":
    pass
