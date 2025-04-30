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
):
    """
    Use the ANARCI program to number the sequence.

    Args:
        seq: An amino acid sequence that you wish to number.

    Returns:
        numbering, chain type, germline information
    """
    from anarci import number as anarci_number

    numbering, chain_type, germline_info = anarci_number(
        seq, allow=allow, assign_germline=True
    )

    if numbering and chain_type in allow:
        return [(_, aa) for _, aa in numbering if aa != "-"], chain_type, germline_info
    elif numbering and chain_type in ["BA", "GD", "AB", "DG"]:
        return (
            [[(_, aa) for _, aa in n if aa != "-"] for n in numbering],
            chain_type,
            germline_info,
        )
    else:
        return False, False, False


def annotate(chain):
    """
    Annotate the sequence of a chain object from TCRDB.TcrPDB
    # e.g. if you have chains B, A and X, you want to force the annotator to return the annotation
    # for B and A but not for X (the antigen)

    returns a dictionary which has the residue ids as key and the annotation as value or is False,
    and chain type which is B/A/G/D/MH1/GA/GB/B2M or False.
    """
    sequence_list, sequence_str = extract_sequence(chain)
    numbering, chain_type, germline_info = call_anarci(sequence_str)

    # Use
    if chain_type:
        chtype = "".join(sorted(chain_type, reverse=True))
    else:
        chtype = False

    if chtype in ("BA", "GD"):
        aligned_numbering = align_scTCR_numbering(
            numbering, sequence_list, sequence_str
        )
        # aligned_numbering = cleanup_scTCR_numbering(aligned_numbering, sequence_list)
        scTCR = True
    # elif chtype == "DC1" or chtype == "RM1":
    #     # Use the scTCR numbering trick; since CD1/MR1 numbering only spans up to residue ~87 and
    #     aligned_numbering = align_scTCR_numbering(
    #         numbering, sequence_list, sequence_str
    #     )
    #     aligned_numbering[0].update(aligned_numbering[1])
    #     aligned_numbering = aligned_numbering[0]  # combine the numbering
    #     aligned_numbering = cleanup_scTCR_numbering(aligned_numbering, sequence_list)
    #     scTCR = False
    else:
        # align the original residue id's to the numbering
        aligned_numbering = align_numbering(numbering, sequence_list)
        scTCR = False

    # aligned numbering is a dictionary of the original residue ids and the new numbering
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


def align_numbering(numbering, sequence_list, alignment_dict={}):
    """
    Align the sequence that has been numbered to the sequence you input.
    The numbered sequence should be "in" the input sequence.
    If not, supply an alignment dictionary.(align sequences and use get_alignment_dict(ali1,ali2))
    """
    if numbering:
        numbered_sequence = "".join([r[1] for r in numbering])
        input_sequence = "".join([r[1] for r in sequence_list])
        if not alignment_dict:
            try:
                numbered_sequence_ali, input_sequence_ali = pairwise_alignment(
                    numbered_sequence, input_sequence
                )
                alignment_dict = get_alignment_dict(
                    input_sequence_ali, numbered_sequence_ali
                )
            except Exception:
                raise Exception(
                    "Could not align numbered sequence to aligned sequence:"
                    + " "
                    + str(numbered_sequence)
                    + " "
                    + str(input_sequence)
                )

        aligned_numbering = {}
        n = -1
        after_flag = False
        for i in range(len(input_sequence)):
            if i in alignment_dict:
                # during
                assert (
                    after_flag is False
                ), "Extra residue in structure than expected from provided sequence"
                assert (
                    input_sequence[i] == numbered_sequence[alignment_dict[i]]
                ), "alignment dictionary failed"
                aligned_numbering[sequence_list[i][0]] = numbering[alignment_dict[i]][0]
                n = numbering[-1][0][0] + 1
            elif n > -1:
                # after
                after_flag = True
                aligned_numbering[sequence_list[i][0]] = (n, " ")
                n += 1
            else:
                # before numbering
                aligned_numbering[sequence_list[i][0]] = ""

        return aligned_numbering
    else:
        return False


def align_scTCR_numbering(numbering, sequence_list, sequence_str):
    """
    Align the sequence that has been numbered to a scTCR structure.
    Args:
        numbering:     numbered list of residues; this is usually a two-element list/tuple from TCRDB.anarci.number
        sequence_list: list of residues (e.g. from a structure) in its original numbering
        sequence_str:  string form of sequence_list
    """
    if numbering:
        numbered_sequence = ["".join([r[1] for r in n]) for n in numbering]
        input_sequence = sequence_str

        aligned_numbering = [{}, {}]

        for ii, a_sequence in enumerate(numbered_sequence):

            # Align each of the joined sequences from the numbering into the target structure sequence in "sequence_str"
            try:
                a_sequence_ali, input_sequence_ali = pairwise_alignment(
                    a_sequence, input_sequence
                )
                alignment_dict = get_alignment_dict(input_sequence_ali, a_sequence_ali)
            except Exception:
                raise Exception(
                    "Could not align numbered sequence to aligned sequence"
                    + "\n"
                    + str(numbered_sequence)
                    + "\n"
                    + str(input_sequence)
                )

            n = -1
            after_flag = False
            # for i in xrange(len(input_sequence)):
            for i in alignment_dict:
                if i in alignment_dict:
                    # during
                    assert (
                        after_flag is False
                    ), "Extra residue in structure than expected from provided sequence"
                    assert (
                        input_sequence[i] == numbered_sequence[ii][alignment_dict[i]]
                    ), "alignment dictionary failed"
                    aligned_numbering[ii][sequence_list[i][0]] = numbering[ii][
                        alignment_dict[i]
                    ][0]
                    n = numbering[ii][-1][0][0] + 1
                elif n > -1:
                    # after
                    after_flag = True
                    aligned_numbering[ii][sequence_list[i][0]] = (n, " ")
                    n += 1
                else:
                    # before numbering
                    aligned_numbering[ii][sequence_list[i][0]] = ""

        return aligned_numbering
    else:
        return False


def cleanup_scTCR_numbering(numbering_dict, sequence_list):
    """
    The scTCR numbering method, while useful for sequences with two domains,
    can have gaps in between (e.g. CD1 molecule of 4lhu).
    This is to close the gaps in the numbering so that residues that were unnumbered by anarci don't move around
    during structural parsing (when they're probably just connections between domains).

    Args:
        numbering_dict: numbered dictionary from align_scTCR_numbering
        sequence_list : sequence list from the structure for alignment.
    """
    positions = [p[0] for p in sequence_list]

    # This gets the last numbered residue in numbering_dict
    lastkey = max(numbering_dict)
    lastidx = positions.index(lastkey)  # Where is this on sequence_list?

    for index in range(1, len(positions)):

        # If we got to the last key, don't bother.
        if index > lastidx:
            break

        key = positions[index]

        # If a target key is not in the numbering dict, see where it fits, then fit a number in it.
        if key not in numbering_dict:

            # Get the left and right bounds of the gap
            left, right = False, False
            lidx, ridx = 0, 0
            lval = (0, " ")
            j = 0

            # Continue iterating left from the missing key until we find one that exists
            while not left:
                key_left = positions[index - j]
                if key_left in numbering_dict:
                    left = True
                    lidx = (
                        index - j
                    )  # Last known index of sequence_list where we know a key exists
                    lval = numbering_dict[key_left]
                else:
                    j += 1

            j = 0
            while not right:
                key_right = positions[index + j]
                if key_right in numbering_dict:
                    right = True
                    ridx = (
                        index + j
                    )  # Last known index of sequence_list on the right where we know a key exists
                else:
                    j += 1

            # For every key between the left and right, fill in
            for k, missing_key in enumerate(positions[lidx + 1 : ridx]):
                numbering_dict[missing_key] = (lval[0] + k + 1, " ")

    return numbering_dict


def get_alignment_dict(ali1, ali2):
    """
    Get a dictionary which tells you the index in sequence 2 that should align with the index in sequence 1 (key)

    ali1:  ----bcde-f---        seq1: bcdef
    ali2:  ---abcd--f---        seq2: abcdf

    alignment_dict={
        0:1,
        1:2,
        2:3,
        4:4
        }

    If the index is aligned with a gap do not include in the dictionary.
    e.g  1 in alignment_dict  --> True
    e.g  3 in alignment_dict  --> False
    """
    assert len(ali1) == len(
        ali2
    ), "aligned sequences must be same lengths (including gaps)"
    alignment_dict = {}
    p1 = -1
    p2 = -1
    for ap in range(len(ali1)):
        if ali1[ap] != "-" and ali2[ap] != "-":
            p1 += 1
            p2 += 1
            alignment_dict[p1] = p2
        elif ali1[ap] != "-":
            p1 += 1
        elif ali2[ap] != "-":
            p2 += 1
    return alignment_dict


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
