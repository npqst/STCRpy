"""
Wrapper for contact_fcc

Calculates residue contacts between different chains
of a single structure.

João Rodrigues @ Utrecht, 2012
"""

import argparse
import os
import sys
from subprocess import Popen, PIPE


def _calculate_contacts(executable: str, pdbfile: str, cutoff: float) -> None:
    """
    Outputs a list of residue contacts based on distance analysis of the PDB file.

    Args:
        executable: path to contact calculation program
        pdbfile: path to PDB-formatted file (.pdb extension)
        cutoff: minimal distance in A to consider a contact between two atoms (float)

    Returns:

    """
    pdbname = os.path.basename(pdbfile)[:-4]
    outfile = os.path.join(os.path.dirname(pdbfile), f"{pdbname}.contacts")
    p = Popen([executable, pdbfile, str(cutoff)], stdout=PIPE, universal_newlines=True)
    p_output, p_error = p.communicate()
    contacts = sorted(list(set([line for line in p_output.split("\n")][:-1])))

    f_contact = open(outfile, "w")
    f_contact.write("\n".join(contacts))
    f_contact.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "files", type=str, nargs="+", help="Pdb file or files to process (or pattern)"
    )
    parser.add_argument(
        "-c", "--cutoff", type=float, default=5.0,
        help="Distance cutoff to evaluate contacts. [Default: 5.0A]"
    )
    parser.add_argument(
        '-e', '--exec', dest="executable", action="store", type=str,
        default="./tools/contact_fcc",
        help="Path to the executable C++ program to calculate contacts [default: ./contact_fcc]"
    )
    args = parser.parse_args()

    # Convert to full paths
    inputs = list(map(os.path.realpath, args.files))

    executable = os.path.abspath(args.executable)
    if not os.path.exists(executable):
        print(f"Path not found: {executable}")
        sys.exit(1)

    for struct in inputs:
        _calculate_contacts(executable, struct, args.cutoff)


if __name__ == "__main__":
    main()
