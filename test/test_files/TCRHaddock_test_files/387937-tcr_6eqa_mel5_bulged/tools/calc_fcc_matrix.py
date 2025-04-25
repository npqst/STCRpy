import argparse
import os
import sys
from time import time, ctime


def parse_contact_file(f_list, ignore_chain) -> list:
    """
    Parses a list of contact files.

    Args:
        f_list:
        ignore_chain:

    Returns:

    """
    if ignore_chain:
        contacts = [[int(line[0:5]+line[6:-1]) for line in open(f)] for f in f_list if f.strip()]
    else:
        contacts = [set([int(line) for line in open(f)]) for f in f_list if f.strip()]
    return contacts


def calculate_fcc(list_a, list_b):
    """
    Calculates the fraction of common elements between two lists taking into account chain IDs.

    Args:
        list_a:
        list_b:

    Returns:

    """
    cc = len(list_a.intersection(list_b))
    cc_v = len(list_b.intersection(list_a))
    return cc, cc_v


def calculate_fcc_nc(list_a, list_b):
    """
    Calculates the fraction of common elements between two lists not taking into account chain IDs.
    Much Slower.

    Args:
        list_a:
        list_b:

    Returns:

    """
    largest, smallest = sorted([list_a, list_b], key=len)
    ncommon = len([ele for ele in largest if ele in smallest])
    return ncommon, ncommon


def calculate_pairwise_matrix(contacts, ignore_chain: bool):
    """
    Calculates a matrix of pairwise fraction of common contacts (FCC).
        Outputs numeric indexes.

        contacts: list_of_unique_pairs_of_residues [set/list]

        Returns pairwise matrix as an iterator, each entry in the form:
        FCC(cplx_1/cplx_2) FCC(cplx_2/cplx_1)

    Args:
        contacts: list_of_unique_pairs_of_residues [set/list]
        ignore_chain:

    Returns: pairwise matrix as an iterator, each entry in the form:
        FCC(cplx_1/cplx_2) FCC(cplx_2/cplx_1).
    """
    contact_lengths = []
    for c in contacts:
        try:
            ic = 1.0/len(c)
        except ZeroDivisionError:
            ic = 0
        contact_lengths.append(ic)

    if ignore_chain:
        calc_fcc = calculate_fcc_nc
    else:
        calc_fcc = calculate_fcc

    for i in range(len(contacts)):
        for k in range(i+1, len(contacts)):
            cc, cc_v = calc_fcc(contacts[i], contacts[k])
            fcc, fcc_v = cc*contact_lengths[i], cc*contact_lengths[k]
            yield i+1, k+1, fcc, fcc_v


def _output_fcc(output, values, f_buffer) -> None:
    """

    Args:
        output:
        values:
        f_buffer:

    Returns:

    """
    buf = []
    for i in values:
        buf.append(i)
        if len(buf) == f_buffer:
            # FIXME: Refactoring
            output("".join([f"{i[0]} {i[1]} {i[2]:1.3f} {i[3]:1.3f}\n" for i in buf]))
            buf = []
    output("".join([f"{i[0]} {i[1]} {i[2]:1.3f} {i[3]:1.3f}\n" for i in buf]))


def main() -> None:
    parser = argparse.ArgumentParser(description="""
Authors: TRELLET Mikael
         RODRIGUES Joao
         MELQUIOND Adrien

         + Undocumented feature: empty contact lists are treated as 0 contacts. (230412 JR)

Calculates a matrix of fraction of common contacts between two or more structures.

WARNING!!
Structures are named numerically in the matrix by the order in which they are read.
By default we sort them by the last number (see _fit*pdb files in it1/analysis/ for example).
""")
    parser.add_argument(
        "contact_files", type=str, nargs="+", help="Contact file or files to process (or pattern)"
    )
    parser.add_argument(
        "-o", "--output", dest="output_file", action="store", type=str, default=sys.stdout,
        help="Output File [default: STDOUT]."
    )
    parser.add_argument(
        "-b", "--buffer_size", dest="buffer_size", action="store", type=int, default=50000,
        help="Number of lines to cache before writing to output file [default: 50000]"
    )
    parser.add_argument(
        "-i", "--ignore_chain", dest="ignore_chain_char", action="store_true",
        help="Ignore chain character in residue code. Use for homomeric complexes."
    )
    parser.add_argument(
        "-H",
        "--haddock-run",
        dest="haddock_run",
        action="store_true",
        help="Touches a final MTX_DONE file. Only relevant for internals of the HADDOCK software."
    )
    args = parser.parse_args()

    # Assume _fit files and sort numerically
    inputs = sorted(args.contact_files, key=lambda x: int(x.split("_")[-1].split(".")[0]))
    if len(inputs) < 2:
        sys.stderr.write("- Provide (at least) two structures to calculate a matrix.")

    sys.stderr.write(f"+ BEGIN: {ctime()}\n")
    if args.ignore_chain_char:
        sys.stderr.write("+ Ignoring chains. Expect a considerable slowdown!!\n")
        exclude_chains = True
    else:
        exclude_chains = False

    t_init = time()
    sys.stderr.write(f"+ Parsing {len(inputs)} contact files\n")

    c = parse_contact_file(inputs, exclude_chains)

    m = calculate_pairwise_matrix(c, exclude_chains)

    if isinstance(args.output_file, str):
        f = open(args.output_file, "w")
    else:
        f = args.output_file

    sys.stderr.write("+ Calculating Matrix\n")
    sys.stderr.write(f"+ Writing matrix to {f.name}\n")
    _output_fcc(f.write, m, args.buffer_size)

    if isinstance(args.output_file, str):
        f.close()
    t_elapsed = time() - t_init
    sys.stderr.write(f"+ END: {ctime()} [{t_elapsed:6.2f} seconds elapsed]\n")
    if args.haddock_run:
        rundir = os.path.dirname(args.output_file)
        touch_file = open(os.path.join(rundir, "MTX_DONE"), "w")
        touch_file.close()


if __name__ == "__main__":
    main()
