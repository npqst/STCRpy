import argparse
import copy
import os
import sys
from time import ctime
from typing import Union


def plugin_core(paramdict, inputlist):
    """

    Args:
        paramdict:
        inputlist:

    Returns:

    """
    print("--> Starting PDB traceback process")

    if inputlist is None:
        traceback = StructureTraceback()
        traceback.setup_run_dir(paramdict)
        traceback.get_start_struc()
        traceback.get_wat_structures()
        traceback.get_it1_structures()
        traceback.get_it0_structures()
        traceback.write_file(verbose=paramdict["verbose"], longout=paramdict["longout"])
    else:
        pdb = []
        for n in inputlist:
            base, extension = os.path.splitext(os.path.basename(n))
            if extension == ".pdb":
                pdb.append(base)
            else:
                try:
                    lines = open(n).readlines()
                    for line in lines:
                        if line == "\n":
                            pass
                        else:
                            pdb.append(line)
                    print(
                        f"    * Supply traceback information for file {os.path.basename(n)}"
                    )
                except Exception:  # TODO: should be more specific
                    raise IOError("    * ERROR: Could not parse file")

        traceback = StructureTraceback()
        traceback.id_structures(filelist=pdb)
        traceback.setup_run_dir(paramdict)
        traceback.get_start_struc()
        traceback.get_wat_structures()
        traceback.get_it1_structures()
        traceback.get_it0_structures()
        traceback.report_query(verbose=paramdict["verbose"])


class StructureTraceback:
    """
    Traceback any structure within a run directory all the way back
    to the individual components that were used in the docking
    """
    def __init__(self) -> None:
        self.fileA_list = []
        self.fileB_list = []
        self.fileC_list = []
        self.fileD_list = []
        self.fileE_list = []
        self.fileF_list = []
        self.complex_list = []
        self.fileit0_list = []
        self.fileit1_list = []
        self.filew_list = []
        self.query = {}
        self.rundir = None
        self.nrstrucit0 = 0
        self.nrstrucit1 = 0
        self.nrstrucw = 0

    @staticmethod
    def _format_line(line: str) -> float:
        """

        Args:
            line:

        Returns:

        """
        base = os.path.splitext(line.strip('""'))
        base2 = (base[0].split(":"))[1]
        base3 = base2.split("_")[-1]

        try:
            return float(base3)
        except Exception:  # TODO: should be more specific
            return float((base3.split("w"))[0])

    @staticmethod
    def _sort_list(inlist: Union[list, None] = None, sortid=None) -> list:
        """

        Args:
            inlist:
            sortid:

        Returns:

        """
        keydir = {}
        for n in inlist:
            keydir[n[sortid]] = n

        name_keys = list(keydir.keys())
        name_keys.sort()

        new = []
        for n in name_keys:
            new.append(keydir[n])

        return new

    @staticmethod
    def _sort_on_index(inlist: Union[list, None] = None, indexlist: Union[list, None] = None,
                       index: bool = False) -> list:
        """

        Args:
            inlist:
            indexlist:
            index:

        Returns:

        """
        indexnr = []
        for n in indexlist:
            indexnr.append(n[0])

        tmp1 = []
        tmplist = copy.deepcopy(inlist)

        if not index:
            for n in indexnr:
                for k in tmplist:
                    if n == k[0]:
                        tmp1.append(k)
                        tmplist.remove(k)

        if index:
            for n in indexnr:
                tmp1.append(inlist[int(n)-1])
            for n in tmplist:
                if n in tmp1:
                    tmplist.remove(n)
                else:
                    pass

        return tmp1+tmplist

    def _fill_list(self) -> None:
        """

        Returns:

        """
        lenit0 = len(self.fileit0_list)
        lenit1 = len(self.fileit1_list)
        lenw = len(self.filew_list)

        tmpit1 = []
        while lenit1 < lenit0:
            tmpit1.append([0.0, 0.0])
            lenit1 = lenit1+1

        tmpw = []
        while lenw < lenit0:
            tmpw.append([0.0, 0.0])
            lenw = lenw+1

        self.fileit1_list = self.fileit1_list+tmpit1
        self.filew_list = self.filew_list+tmpw

    def id_structures(self, filelist: Union[list, None] = None) -> None:
        """

        Args:
            filelist:

        Returns:

        """
        currdir = os.getcwd()
        base, ext = os.path.split(currdir)

        if ext == "it0":
            lib = "it0"
        elif ext == "it1":
            lib = "it1"
        elif ext == "water":
            lib = "water"
        else:
            print("    * ERROR: Structure not present in either it0, it1 or water directory")
            sys.exit(0)

        tmp = []
        for n in filelist:
            try:
                base = n.split("_")[1]
                # TODO: check on example as previous version "float(base)[0]" seems to be faulty
                tmp.append(float(base[0]))
            except Exception:  # TODO: should be more specific
                base = n.split("_")[1]
                tmp.append(float(base.split("w")[0]))

        self.query[lib] = tmp

    def setup_run_dir(self, paramdict: dict) -> None:
        """
        Setup the running directory from user input (inputdir) or as
        the current working directory if no input.

        Args:
            paramdict:

        Returns:

        """
        if paramdict["inputdir"] and os.path.exists(paramdict["inputdir"]):
            if not os.path.exists(os.path.join(os.path.abspath(paramdict["inputdir"]), "begin")):
                print(f"    * ERROR: No begin directory in {paramdict['inputdir']}")
                sys.exit(0)
            else:
                self.rundir = os.path.abspath(paramdict["inputdir"])
        elif not paramdict["inputdir"]:
            self.rundir = os.getcwd()
        else:
            print(f"    * ERROR: Could not find directory {paramdict['inputdir']}")
            sys.exit(0)
        print(f"    * Working directory: {self.rundir}")

    def get_start_struc(self) -> None:
        """
        Getting all starting structures from the file_X.list files in the begin directory.

        Returns:

        """
        begindir = os.path.join(self.rundir, "begin")
        os.chdir(begindir)

        files = [
            "file_A.list", "file_B.list", "file_C.list", "file_D.list", "file_E.list", "file_F.list"
        ]
        files2 = [
            "fileA_list", "fileB_list", "fileC_list", "fileD_list", "fileE_list", "fileF_list"
        ]

        for file_list in files:
            if os.path.isfile(file_list):
                lines = open(file_list).readlines()

                for line in lines:
                    tmp = (line.split("/"))[-1].strip('"\n')
                    # TODO: check how tmp is used below as in the previous version
                    # the variable hasn't been used as it was kept as string literal
                    exec(f"self.{files2[files.index(file_list)]}.append({tmp})")

        """"Make combinations equal to generate_complex.inp and store in complex_list"""
        # TODO: this needs to be simplified
        if len(self.fileA_list) > 0:
            for structureA in self.fileA_list:
                if len(self.fileB_list) > 0:
                    for structureB in self.fileB_list:
                        if len(self.fileC_list) > 0:
                            for structureC in self.fileC_list:
                                if len(self.fileD_list) > 0:
                                    for structureD in self.fileD_list:
                                        if len(self.fileE_list) > 0:
                                            for structureE in self.fileE_list:
                                                if len(self.fileF_list) > 0:
                                                    for structureF in self.fileF_list:
                                                        self.complex_list.append(
                                                            structureA + "/" + structureB + "/" +
                                                            structureC + "/" + structureD+"/" +
                                                            structureE + "/" + structureF
                                                        )
                                                else:
                                                    self.complex_list.append(
                                                        structureA + "/" + structureB + "/" +
                                                        structureC + "/" + structureD + "/" +
                                                        structureE
                                                    )
                                        else:
                                            self.complex_list.append(
                                                structureA + "/" + structureB + "/" + structureC +
                                                "/" + structureD
                                            )
                                else:
                                    self.complex_list.append(
                                        structureA + "/" + structureB + "/" + structureC
                                    )
                        else:
                            self.complex_list.append(structureA + "/" + structureB)
                else:
                    self.complex_list.append(structureA)
        else:
            print("    * ERROR No starting structures found in the begin directory")
            sys.exit(0)

    def get_it0_structures(self) -> None:
        """

        Returns:

        """
        begindir = os.path.join(self.rundir, "structures", "it0")
        os.chdir(begindir)

        if os.path.isfile("file.list"):
            lines = open("file.list").readlines()

            for line in lines:
                if line == '\n':
                    pass
                else:
                    tmp = [self._format_line(line.split()[0]), float(line.split()[2])]
                    self.fileit0_list.append(tmp)
        else:
            print("    * ERROR: No file.list found in it0 directory. Nothing to trace means stop")
            sys.exit(0)

        self.nrstrucit0 = len(self.fileit0_list)
        nrstrucbg = len(self.complex_list)

        self.complex_list = self.complex_list * int((round(self.nrstrucit0/nrstrucbg))+1)
        self.complex_list = self.complex_list[0:self.nrstrucit0]

        # first sort on structure number low->high
        self.fileit0_list = self._sort_list(inlist=self.fileit0_list, sortid=0)

        for n in range(len(self.complex_list)):  # Match complex_list to it0 structures
            self.fileit0_list[n].append(self.complex_list[n])

        # then sort on HADDOCK score low->high
        self.fileit0_list = self._sort_list(inlist=self.fileit0_list, sortid=1)
        # Sort according to index it1
        self.fileit0_list = self._sort_on_index(
            inlist=self.fileit0_list, indexlist=self.fileit1_list, index=True
        )

    def get_it1_structures(self) -> None:
        """

        Returns:

        """
        begindir = os.path.join(self.rundir, "structures", "it1")
        os.chdir(begindir)

        if os.path.isfile("file.list"):
            lines = open("file.list").readlines()

            for line in lines:
                if line == "\n":
                    pass
                else:
                    tmp = [self._format_line(line.split()[0]), float(line.split()[2])]
                    self.fileit1_list.append(tmp)
        else:
            print("    * ERROR: No file.list found in it1 directory. Nothing to trace means stop")
            sys.exit(0)

        self.nrstrucit1 = len(self.fileit1_list)
        # Sort according to index W. This matches it1 to W
        self.fileit1_list = self._sort_on_index(
            inlist=self.fileit1_list, indexlist=self.filew_list, index=False
        )

    def get_wat_structures(self) -> None:
        """

        Returns:

        """
        begindir = os.path.join(self.rundir, "structures", "it1", "water")
        os.chdir(begindir)

        if os.path.isfile("file.list_all"):
            lines = open("file.list_all").readlines()
        else:
            lines = open("file.list").readlines()

        for line in lines:
            if line == "\n":
                pass
            else:
                tmp = [self._format_line(line.split()[0]), float(line.split()[2])]
                self.filew_list.append(tmp)

        self.nrstrucw = len(self.filew_list)

        if self.nrstrucw > 0:
            # Sort on HADDOCK score low->high.
            self.filew_list = self._sort_list(inlist=self.filew_list, sortid=1)
        else:
            print("    No file.list of file.list_all found in water refinement directory. \
            Only traceback from it1 to it0")

    def write_file(self, verbose: bool = False, longout: bool = False) -> None:
        """

        Args:
            verbose:
            longout:

        Returns:

        """
        os.chdir(self.rundir)

        if verbose:
            outfile = sys.stdout
        else:
            outfile = open("traceback.list", "w")
            print(
                (
                    "    * Traceback information written to file 'traceback.list'"
                    f" in directory {self.rundir}")
            )

        outfile.write("*" * 113 + "\n")
        outfile.write(f"Structure traceback information for run {self.rundir}\n")
        outfile.write(f"Date/time: {ctime()}\n")
        outfile.write(
            "Number of structures: {self.nrstrucit0} in it0, {self.nrstrucit1} in it1 and {self.nrstrucw} in water refinement\n"
        )
        outfile.write(
            (
                "Sorting order: water(struct. nr.) "
                "matches it1 (struct. nr.) "
                "matches it0 (struct. nr.) matches         input structures.\n"
            )
        )
        outfile.write("*" * 113 + "\n")
        outfile.write(
            (
                "      complex                             "
                "it0      hscoreit0       it1      hscoreit1"
                "              water      hscorew\n"
            )
        )

        # TODO: remove code duplication below with a function
        if longout:
            self._fill_list()
            for n in range(len(self.fileit0_list)):
                # FIXME: refactor
                outfile.write(
                    "{:>35}{:10.0f}{:15.4f}{:10.0f}{:15.4f}{:10.0f}{:15.4f}\n".format(
                        self.fileit0_list[n][2],
                        self.fileit0_list[n][0],
                        self.fileit0_list[n][1],
                        self.fileit1_list[n][0],
                        self.fileit1_list[n][1],
                        self.filew_list[n][0],
                        self.filew_list[n][1]
                    )
                )
        else:
            for n in range(len(self.filew_list)):
                # FIXME: refactor
                outfile.write(
                    "{:>35}{:10.0f}{:15.4f}{:10.0f}{:15.4f}{:10.0f}{:15.4f}\n".format(
                        self.fileit0_list[n][2],
                        self.fileit0_list[n][0],
                        self.fileit0_list[n][1],
                        self.fileit1_list[n][0],
                        self.fileit1_list[n][1],
                        self.filew_list[n][0],
                        self.filew_list[n][1]
                    )
                )

        if not verbose:
            outfile.close()

    def report_query(self, verbose: bool = False) -> None:
        """

        Args:
            verbose:

        Returns:

        """
        os.chdir(self.rundir)

        if verbose:
            outfile = sys.stdout
        else:
            outfile = open("traceback.list", "w")
            print(
                f"    * Traceback information written to file 'traceback.list' in directory {self.rundir}"  # noqa
            )

        outfile.write("*" * 113 + "\n")
        outfile.write(f"Structure traceback information for run {self.rundir}\n")
        outfile.write(f"Date/time: {ctime()}\n")
        outfile.write(
            "Number of structures: {:d} in it0, {:d} in it1 and {:d} in water refinement\n".format(
                self.nrstrucit0, self.nrstrucit1, self.nrstrucw
            )
        )
        outfile.write(
            (
                "Sorting order: water(struct. nr.) "
                "matches it1 (struct. nr.) "
                "matches it0 (struct. nr.) matches         input structures.\n"
            )
        )

        outfile.write("*" * 113 + "\n")
        outfile.write(
            (
                "      complex                             "
                "it0      hscoreit0       it1      hscoreit1"
                "              water      hscorew\n"
            )
        )

        if self.query.keys() == ["water"]:
            for n in self.query["water"]:
                for k in self.filew_list:
                    if n == k[0]:
                        outfile.write(
                            "{:>35}{:10.0f}{:15.4f}{:10.0f}{:15.4f}{:10.0f}{:15.4f}\n".format(
                                self.fileit0_list[self.filew_list.index(k)][2],
                                self.fileit0_list[self.filew_list.index(k)][0],
                                self.fileit0_list[self.filew_list.index(k)][1],
                                self.fileit1_list[self.filew_list.index(k)][0],
                                self.fileit1_list[self.filew_list.index(k)][1],
                                n, k[1]
                            )
                        )
        elif self.query.keys() == ['it1']:
            for n in self.query['it1']:
                for k in self.fileit1_list:
                    if n == k[0]:
                        outfile.write(
                            "{:>35}{:10.0f}{:15.4f}{:10.0f}{:15.4f}{:10.0f}{:15.4f}\n".format(
                                self.fileit0_list[self.fileit1_list.index(k)][2],
                                self.fileit0_list[self.fileit1_list.index(k)][0],
                                self.fileit0_list[self.fileit1_list.index(k)][1], n, k[1],
                                self.filew_list[self.fileit1_list.index(k)][0],
                                self.filew_list[self.fileit1_list.index(k)][1]
                            )
                        )
        elif self.query.keys() == ["it0"]:
            self._fill_list()
            for n in self.query["it0"]:
                for k in self.fileit0_list:
                    if n == k[0]:
                        outfile.write(
                            "{:>35}{:10.0f}{:15.4f}{:10.0f}{:15.4f}{:10.0f}{:15.4f}\n".format(
                                k[2], n, k[1], self.fileit1_list[self.fileit0_list.index(k)][0],
                                self.fileit1_list[self.fileit0_list.index(k)][1],
                                self.filew_list[self.fileit0_list.index(k)][0],
                                self.filew_list[self.fileit0_list.index(k)][1]
                            )
                        )

        if not verbose:
            outfile.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="""
    ==========================================================================================

    Author:			Marc van Dijk, Department of NMR spectroscopy, Bijvoet Center
                for Biomolecular Research, Utrecht university, The Netherlands
    Copyright (C):		2006 (DART project)
    DART version:		1.0  (01-01-2007)
    DART plugin: 		PDBtraceback.py
    Input:			Nothing, PDB files or file-lists and any of the allowed options
                any order and combined
    Output:			A file called traceback.list or terminal output
    Plugin excecution:	Either command line driven (use -h/--help for the option) or as
                part of a DART batch sequence. Excecute anywhere inside the run
                directory of choice.
    Plugin function:	This script allows you to "traceback" any structure anywhere
                within a run directory all the way back to the structures of the
                individual components used in the docking of that particular
                structure.
    Examples:		PDBtraceback.py
                PDBtraceback.py -f test.pdb
    Plugin dependencies:	None

    for further information, please contact:
                - DART website (http://www.nmr.chem.uu.nl/DART)
                - email: abonvin@chem.uu.nl

    If you are using this software for academic purposes please quoting the following
    reference:

    ==========================================================================================
    """)

    # TODO: check if argument parsing works as expected after move to argparse for various use cases
    parser.add_argument("inputs", type=str, nargs="+", help="Supply pdb or file.nam inputfile(s). \
    Standard UNIX selection syntax accepted")
    parser.add_argument(
        "-d", "--dir", dest="inputdir", type=str, help="Directory path of HADDOCK run."
    )
    parser.add_argument(
        "-l", "--longout", action="store_true", dest="longout", default=False,
        help="Long output also includes all structures rejected in it0 and it1, default=False."
    )
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", default=False,
                        help="All output to standard output.")
    args = parser.parse_args()

    params_dict = dict()
    params_dict["inputdir"] = args.inputdir
    params_dict["input"] = [os.path.realpath(i) for i in args.inputs]
    params_dict["longout"] = args.longout
    params_dict["verbose"] = args.verbose

    """Envoce main functions"""
    plugin_core(params_dict, inputlist=params_dict['input'])

    """Say goodbye"""
    print("--> Thanks for using PDBtraceback, bye")


if __name__ == "__main__":
    main()
