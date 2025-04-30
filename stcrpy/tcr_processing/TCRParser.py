"""
Created on 3 April 2024
@author: Nele Quast, based on leem

TCRParser object which is based on ABDB's AntibodyParser and BioPython's PDB parser.
"""

from itertools import combinations, product
import sys
import os
from collections import defaultdict
import warnings

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import NeighborSearch

# TCRDB
from .annotate import annotate, extract_sequence, align_numbering

from ..utils.error_stream import ErrorStream

from .TCRStructure import TCRStructure
from .Model import Model
from .TCR import TCR, abTCR, gdTCR
from .MHC import MHC, MH1, MH2, CD1, MR1, scMH1, scCD1, scMH2
from .Holder import Holder
from .TCRchain import TCRchain
from .MHCchain import MHCchain
from .AGchain import AGchain
from Bio.PDB.Residue import Residue
from .Fragment import Fragment
from .Chemical_components import is_aa, is_common_buffer, get_res_type, is_carbohydrate

MHC_CUTOFF = {
    "MH1:B2M": [32, 37, 32, 32],
    "CD1:B2M": [32, 37, 32, 32],
    "GA1L:B2M": [33, 37, 33, 32],           # added for GA1L case 2po6
    "MR1:B2M": [32, 37, 32, 32],
    "GA1:B2M": [32, 37, 32, 32],
    "GB:GA": [22, 32, 35, 29],
}


class TCRParser(PDBParser, MMCIFParser):
    def __init__(self, PERMISSIVE=True, get_header=True, QUIET=False):
        """
        Initialise the PDB parser. This is currently set to using IMGT's numbering scheme and uses the IMGT-defined CDRs.

        """
        self.pdb_parser = PDBParser(PERMISSIVE, get_header, None, QUIET)
        self.mmcif_parser = MMCIFParser(None, QUIET)
        self.QUIET = QUIET
        # Structures are numbered using anarci.
        self.numbering_method = "anarci"

        # Choose the numbering scheme and CDR definition, though by default we'll use the IMGT schemes. Have these down for reference use.
        self.numbering_scheme = "imgt"
        self.definition = "imgt"

    def _create_chain(self, chain, new_chain_id, numbering, chain_type):
        """
        Create a new TCR or MHC chain.
        Residues before the numbered region are now ignored.
        """
        if chain_type in ["D", "A", "B", "G"]:
            newchain = TCRchain(new_chain_id)
        elif chain_type in [
            "MH1",
            "CD1",
            "B2M",
            "GA",
            "GB",
            "GA1L",
            "GA2L",
            "GA1",
            "GA2",
            "MR1",
        ]:
            newchain = MHCchain(new_chain_id)

        newchain.numbering = numbering
        unnumbered_list = []
        added = False

        for residue in chain.get_list():
            # check whether to add the residue to the new chain (have we got numbering for it)
            add = False
            if residue.id in numbering:
                if numbering[residue.id]:
                    add = True
                    res_id = (
                        residue.id[0],
                        numbering[residue.id][0],
                        numbering[residue.id][1],
                    )  # field and reannotated things.

            # if we should add it, add it
            if add:
                added = True
                newresidue = Residue(res_id, residue.resname, residue.segid)
                for atom in residue.get_list():
                    newresidue.add(atom.copy())
                newresidue.imgt_numbered = True
                newchain.add(newresidue)

            # else add it to the unnumbered list - this will include the HETATOMs - analyse them to find haptens.
            elif added:
                unnumbered_list.append(residue)

        # add the unnumbered residues into the chain - renumbering so that they follow on from the numbered regions.
        ended = sorted([i for i in numbering.values() if i != ""])[-1][0]
        for residue in unnumbered_list:
            ended += 1
            res_id = (residue.id[0], ended, " ")
            newresidue = Residue(res_id, residue.resname, residue.segid)
            for atom in residue.get_list():
                newresidue.add(atom.copy())
            newchain.add(newresidue)
            newchain.add_unnumbered(newresidue)

        newchain.analyse(chain_type)
        return newchain

    def _create_scTCR_chains(
        self, chain, new_chain_id, numbering_1, numbering_2, chain_type1, chain_type2
    ):
        """
        Create two TCR chains to be paired up as a TCR.
        This is effectively a TCR chain with a modified unnumbered list generation
        Residues before/after the numbered region are now ignored.
        """
        newchain1 = TCRchain(new_chain_id.lower())
        newchain2 = TCRchain(new_chain_id.upper())

        if chain_type2 in ["G", "A"]:
            # Just a trick to
            return self._create_scTCR_chains(
                chain, new_chain_id, numbering_2, numbering_1, chain_type2, chain_type1
            )

        newchain1.numbering = numbering_1
        newchain2.numbering = numbering_2
        newchains = [newchain1, newchain2]

        unnumbered_set = set()
        added = False
        numbered_pos = set(numbering_1.keys()) | set(numbering_2.keys())

        for i, numbering in enumerate([numbering_1, numbering_2]):
            # Get the ith newchain
            newchain = newchains[i]

            for residue in chain.get_list():
                # check whether to add the residue to the new chain (have we got numbering for it)
                add = False
                if residue.id in numbering:
                    if numbering[residue.id]:
                        add = True
                        res_id = (
                            residue.id[0],
                            numbering[residue.id][0],
                            numbering[residue.id][1],
                        )  # field and reannotated things.

                # if we should add it, add it
                if add:
                    added = True
                    newresidue = Residue(res_id, residue.resname, residue.segid)
                    for atom in residue.get_list():
                        newresidue.add(atom.copy())
                    newresidue.imgt_numbered = True
                    newchain.add(newresidue)

                # else add it to the unnumbered list - this will include the HETATOMs - analyse them to find haptens.
                elif added and residue.id not in numbered_pos:
                    unnumbered_set.add(residue)

        # add the unnumbered residues into the chain - renumbering so that they follow on from the numbered regions.
        ended = sorted(numbering_1.values())[-1][0]  # get the last numbered value.

        for residue in sorted(unnumbered_set, key=lambda z: z.id[1]):
            ended += 1
            res_id = (residue.id[0], ended, " ")
            newresidue = Residue(res_id, residue.resname, residue.segid)
            for atom in residue.get_list():
                newresidue.add(atom.copy())

            newchain1.add(newresidue)
            newchain1.add_unnumbered(newresidue)

        newchain1.analyse(chain_type1)
        newchain2.analyse(chain_type2)

        return newchain1, newchain2

    def get_tcr_structure(
        self, id, file, prenumbering=None, ali_dict={}, crystal_contacts=[]
    ):
        """
        Post processing of the TCRPDB.Bio.PDB structure object into a TCR context.

        id: a string to identify the structure
        file: the path to the .pdb file

        optional:
            prenumbering: prenumbering for the chains in the structure.
        """
        self.warnings = ErrorStream()
        # get a structure object from biopython.
        _, ext = os.path.splitext(file)
        if ext.lower() == ".pdb":
            structure = self.pdb_parser.get_structure(id, file)
            self.current_parser = self.pdb_parser
        elif ext.lower() in [".cif", ".mmcif"]:
            structure = self.mmcif_parser.get_structure(id, file)
            self.current_parser = self.mmcif_parser
        else:
            self.warnings.write(f"Unrecognised structure file format: {file}")
            raise ValueError

        # Create a new TCRStructure object
        tcrstructure = TCRStructure(structure.id)

        # Set and analyse header information
        tcrstructure.set_header(structure.header)
        self._analyse_header(tcrstructure)

        # iterate over the models in the structure
        # iterate backwards through the model list - delete old structure as we go
        # e.g. NMR structures will be extremely memory expensive (72 models!)

        for mid in range(len(structure.child_list) - 1, -1, -1):
            # add a model to the TCR structure
            model = structure.child_list[mid]
            newmodel = Model(model.id)
            tcrstructure.add(newmodel)

            # initialise holder objects for holding TCR, MHC and non-TCR/non-MHC (antigen) chains.
            agchains = Holder("Antigen")
            trchains = Holder("TCRchain")
            mhchains = Holder("MHCchain")
            newmodel.add(agchains)
            newmodel.add(trchains)
            newmodel.add(mhchains)

            # iterate over the chains in the model
            for chain in model.get_list():
                # try to number the sequence found in the structure
                if prenumbering and chain.id in prenumbering:
                    if len(prenumbering[chain.id]) == 2:
                        numbering = [{}, {}]
                        region_types = ["", ""]

                        numbering[0], region_types[0] = self._prenumbered(
                            chain, prenumbering, ali_dict, n=0
                        )
                        numbering[1], region_types[1] = self._prenumbered(
                            chain, prenumbering, ali_dict, n=1
                        )
                        rtypes = sorted(region_types)

                        # Check that we have a beta/alpha domain or gamma/delta domain
                        if rtypes == ["A", "B"] or rtypes == ["D", "G"]:
                            chain_type = "".join(region_types)
                            scTCR = True
                        # if not, just take the first region and warn the user
                        else:
                            chain_type = region_types[0]
                            numbering = numbering[0]
                            scTCR = False
                            print(
                                "Warning multiple variable regions of the same type (%s) found on chain %s.\nTaking the first variable region only."
                                % (chain_type, chain.id),
                                file=self.warnings,
                            )

                    elif prenumbering[chain.id][0][-1] not in ["B", "A", "D", "G"]:
                        numbering, chain_type, scTCR = annotate(chain)

                    else:
                        numbering, chain_type = self._prenumbered(
                            chain, prenumbering, ali_dict, n=0
                        )
                        scTCR = False

                else:
                    numbering, chain_type, germline_info, scTCR = annotate(chain)

                if chain.id in tcrstructure.header["chain_details"]:  # clean this up!!!
                    engineered = tcrstructure.header["chain_details"][chain.id][
                        "engineered"
                    ]
                    details = tcrstructure.header["chain_details"][chain.id]
                else:
                    engineered = False
                    details = {"molecule": "unknown", "engineered": False}

                details["genetic_origin"] = germline_info

                if numbering and chain_type in ["G", "D", "B", "A"]:
                    # create a new TCR chain
                    newchain = self._create_chain(
                        chain, chain.id, numbering, chain_type
                    )
                    newchain.set_engineered(engineered)
                    newchain.xtra.update(details)
                    trchains.add(newchain)

                elif numbering and chain_type in [
                    "MH1",
                    "CD1",
                    "GA",
                    "GB",
                    "B2M",
                    "GA1L",
                    "GA2L",
                    "GA1",
                    "GA2",
                    "MR1",
                ]:
                    newchain = self._create_chain(
                        chain, chain.id, numbering, chain_type
                    )
                    newchain.set_engineered(engineered)
                    newchain.xtra.update(details)
                    mhchains.add(newchain)

                elif numbering and scTCR:
                    # Separate numbering into two domains
                    types = list(chain_type)
                    domain1, domain2 = numbering

                    chain1, chain2 = self._create_scTCR_chains(
                        chain, chain.id, domain1, domain2, types[0], types[1]
                    )
                    chain1.set_engineered(engineered)
                    chain1.xtra.update(details)
                    chain2.set_engineered(engineered)
                    chain2.xtra.update(details)

                    # We know this is a TCR -- except for 2p1y.
                    if (
                        chain1.child_dict[(" ", 104, " ")]["CA"]
                        - chain2.child_dict[(" ", 104, " ")]["CA"]
                    ) <= 22:
                        obs_chaintypes = set([chain1.chain_type, chain2.chain_type])
                        if not obs_chaintypes - set(["A", "B"]):
                            tcr = abTCR(chain1, chain2)
                        elif not obs_chaintypes - set(["G", "D"]):
                            tcr = gdTCR(chain1, chain2)
                        elif not obs_chaintypes - set(["B", "D"]):
                            tcr = abTCR(chain1, chain2)         # initial way to deal with narci missclassification of alpha chains as delta chains
                            # tcr = dbTCR(chain1, chain2)

                        tcr.scTCR = True  #
                        newmodel.add(tcr)
                        if chain1.id in trchains:
                            trchains.detach_child(chain1.id)
                        if chain2.id in trchains:
                            trchains.detach_child(chain2.id)

                    else:
                        trchains.add(chain1)
                        trchains.add(chain2)

                # add chain to "antigen" chains
                else:
                    newchain = self._create_ag_chain(chain)
                    newchain.set_engineered(engineered)
                    newchain.xtra.update(details)
                    agchains.add(newchain)

            # try to pair the TCR chains to form TCRs. Use a heuristic for now.
            if not scTCR:
                pairings = self._pair_chains(trchains)

                for pair in pairings:
                    trchains.detach_child(pair[0].id)
                    trchains.detach_child(pair[1].id)

                    obs_chaintypes = set([pair[0].chain_type, pair[1].chain_type])
                    if not obs_chaintypes - set(["A", "B"]):
                        tcr = abTCR(pair[0], pair[1])
                    elif not obs_chaintypes - set(["G", "D"]):
                        tcr = gdTCR(pair[0], pair[1])
                    elif not obs_chaintypes - set(["B", "D"]):
                        # tcr = dbTCR(pair[0], pair[1])
                        tcr = abTCR(pair[0], pair[1])

                    else:
                        self.warnings.write(
                            "Unusual pairing between %s (V%s) and %s (V%s) has been detected. Treating as separate TCR chains.\n"
                            % (
                                pair[0].id,
                                pair[0].chain_type,
                                pair[1].id,
                                pair[1].chain_type,
                            )
                        )
                        trchains.add(pair[0])
                        trchains.add(pair[1])
                        continue
                    newmodel.add(tcr)

            elif scTCR and trchains:
                pairings = self._pair_chains(trchains)

                for pair in pairings:
                    trchains.detach_child(pair[0].id)
                    trchains.detach_child(pair[1].id)

                    obs_chaintypes = set([pair[0].chain_type, pair[1].chain_type])
                    if not obs_chaintypes - set(["A", "B"]):
                        tcr = abTCR(pair[0], pair[1])
                    elif not obs_chaintypes - set(["G", "D"]):
                        tcr = gdTCR(pair[0], pair[1])
                    elif not obs_chaintypes - set(["B", "D"]):
                        tcr = abTCR(pair[0], pair[1])
                        # tcr = dbTCR(pair[0], pair[1])
                    else:
                        self.warnings.write(
                            "Unusual pairing between %s (V%s) and %s (V%s) has been detected. Treating as separate TCR chains.\n"
                            % (
                                pair[0].id,
                                pair[0].chain_type,
                                pair[1].id,
                                pair[1].chain_type,
                            )
                        )
                        trchains.add(pair[0])
                        trchains.add(pair[1])
                        continue

                    newmodel.add(tcr)

            # Pair up the MHC chains -- whether that's GA and GB or MH1 with B2M
            pairings = self._pair_mhc(mhchains)
            for pair in pairings:
                mhchains.detach_child(pair[0].id)
                mhchains.detach_child(pair[1].id)

                obs_chaintypes = set([pair[0].chain_type, pair[1].chain_type])
                if (
                    not (obs_chaintypes - set(["MH1", "B2M"]))
                    or not (obs_chaintypes - set(["GA1", "GA2"]))
                    or not (obs_chaintypes - set(["GA1", "B2M"]))
                ):
                    mhc = MH1(pair[0], pair[1])
                elif not (obs_chaintypes - set(["GA", "GB"])):
                    mhc = MH2(pair[0], pair[1])
                elif not (obs_chaintypes - set(["CD1", "B2M"])) or not (
                    obs_chaintypes - set(["GA1L", "GA2L"])) or not (
                        obs_chaintypes - set(["GA1L", "B2M"])
                ):
                    mhc = CD1(pair[0], pair[1])
                elif not (obs_chaintypes - set(["MR1", "B2M"])):
                    mhc = MR1(pair[0], pair[1])
                else:
                    raise ValueError(f'MHC pairing {pair} could not be assigned.')

                newmodel.add(mhc)

            # allow instantiation of single chain MH1 type MH class if the alpha helices forming chain has been observed
            # allow instantiation of single chain MH2 type MH class if one of the GA or GB chain has been observed
            ids_to_detach = []
            for mhc_chain in mhchains:
                if mhc_chain.chain_type in ["MH1", "GA1", "GA2"]:
                    ids_to_detach.append(mhc_chain.id)
                    sc_mhc = scMH1(mhc_chain)
                    newmodel.add(sc_mhc)
                elif mhc_chain.chain_type in ["CD1", "GA1L"]:
                    ids_to_detach.append(mhc_chain.id)
                    sc_mhc = scCD1(mhc_chain)
                    newmodel.add(sc_mhc)
                elif mhc_chain.chain_type in ["GA", "GB"]:
                    ids_to_detach.append(mhc_chain.id)
                    sc_mhc = scMH2(mhc_chain)
                    newmodel.add(sc_mhc)
                    warnings.warn(
                        f"Single chain MH class II instantiated with chain type {mhc_chain.chain_type}. It is possible the other MHC class II chain has not been identified."
                    )

            for mhc_chain_id in ids_to_detach:
                mhchains.detach_child(mhc_chain_id)

            # Match MHC+antigen complex with a TCR
            self._match_units(newmodel, trchains, mhchains, agchains, crystal_contacts)
            del structure.child_list[
                mid
            ]  # delete the structure model list (goes backwards so indexing is not affected)

            # Delete empty holders
            empty_holders = [
                holder.id for holder in newmodel.child_list if not holder.child_list
            ]
            for holder_id in empty_holders:
                newmodel.detach_child(holder_id)

        del structure
        if not self.QUIET and self.warnings.log:
            sys.stderr.write("\n".join(self.warnings.log))
            sys.stderr.write("\n")
        tcrstructure.warnings = self.warnings

        return tcrstructure

    def _analyse_header(self, header):
        """
        Analysis of the header that has been parsed by Biopython
        We add information for the various chains and have a look for engineered and hapten flags.
        Add more information to this parser.
        """
        if isinstance(header, TCRStructure):
            header = header.get_header()
        elif not header:
            header = {}

        header["chain_details"] = {}
        if "compound" in header:
            for compound in header["compound"]:
                # iteration over details.
                if "chain" in header["compound"][compound]:
                    # get the chains that the compound is refering to.
                    chains = [
                        c.strip().upper()
                        for c in header["compound"][compound]["chain"].split(",")
                        if len(c.strip()) == 1
                    ]

                    for chain in chains:
                        if chain not in header["chain_details"]:
                            header["chain_details"][chain] = {}

                    if "molecule" in header["compound"][compound]:
                        # add molecule annotation to each chain
                        for chain in chains:
                            header["chain_details"][chain]["molecule"] = header[
                                "compound"
                            ][compound]["molecule"]
                    else:
                        for chain in chains:
                            header["chain_details"][chain]["molecule"] = "unknown"

                    if "engineered" in header["compound"][compound]:
                        if (
                            "no" in header["compound"][compound]["engineered"]
                            or "false" in header["compound"][compound]["engineered"]
                            or not header["compound"][compound]["engineered"]
                        ):
                            header["compound"][compound]["engineered"] = False
                        else:
                            header["compound"][compound]["engineered"] = True
                        for chain in chains:
                            header["chain_details"][chain]["engineered"] = header[
                                "compound"
                            ][compound]["engineered"]
                    else:
                        for chain in chains:
                            header["chain_details"][chain]["engineered"] = False
                else:
                    continue

            # analyse the journal reference and the title for references to hapten or scfv
            # compile title-like text
            title = (
                header["journal_reference"].lower()
                + " ".join(header["structure_reference"]).lower()
            )
            if "journal" in header:
                title += header["journal"].lower()
        else:
            sys.stderr.write("Header could not be parsed")

    def _create_ag_chain(self, chain):
        """
        Create a new 'antigen' chain - this just means it is not a TCR chain.
        """
        newchain = AGchain(chain.id)
        for residue in chain.get_list():
            newresidue = Residue(residue.id, residue.resname, residue.segid)
            newchain.add(newresidue)
            for atom in residue.get_list():
                newresidue.add(atom.copy())
        newchain.set_type()
        return newchain

    def _pair_chains(self, chains):
        """
        Method to pair beta/alpha and gamma/delta chains to form TCRs.
        Currently this is based off of ABDB.AbPDB's chain pairing method where the
        distance between positions 104 are calculated using the same 22A cutoff.
        This is a simple heuristic for now.
        """
        pairings = []
        # We use a known distance between conserved cysteine residues at the interface
        points = {
            "B": (" ", 104, " "),
            "A": (" ", 104, " "),
            "D": (" ", 104, " "),
            "G": (" ", 104, " "),
        }

        for pair in combinations(chains, 2):
            if pair[0].chain_type != pair[1].chain_type:
                try:
                    a1 = pair[0].child_dict[points[pair[0].chain_type]].child_dict["CA"]
                    a2 = pair[1].child_dict[points[pair[1].chain_type]].child_dict["CA"]
                except KeyError:
                    continue
                if a1 - a2 < 22:
                    pairings.append(pair)
        return pairings

    def _pair_mhc(self, chains):
        """
        This is a heuristic that pairs MHC chains together. In theory, we should have a GA-GB chain (MHC2) or a GA1/GA2-C (MHC1).
        Use an arbitrary cutoff of 45A for pairing the MHC for now.
        Where possible, use the conserved cysteine in the GA2/GB domains (resi 1074 or 11) and pair up with another
        conserved point (Cys 104 in B2M, N86 in GA)

        Impose an angle cutoff so that the cysteine of the B2M points in the correct orientation
        """
        pairings = []
        points = {
            "MH1": [(" ", 15, " "), (" ", 51, " ")],
            "GA1": [(" ", 15, " "), (" ", 51, " ")],
            "GA2": [(" ", 15, " "), (" ", 51, " ")],
            "CD1": [(" ", 15, " "), (" ", 51, " ")],  # pretty similar to MH1
            "GA1L": [(" ", 15, " "), (" ", 51, " ")],  # pretty similar to MH1
            "MR1": [(" ", 15, " "), (" ", 51, " ")],  # pretty similar to MH1
            "B2M": [(" ", 23, " "), (" ", 104, " ")],
            "GA": [(" ", 29, " "), (" ", 37, " ")],
            "GB": [(" ", 39, " "), (" ", 64, " ")],
        }

        acceptable_types = [
            "MH1:B2M",
            "GB:GA",
            "CD1:B2M",
            "MR1:B2M",
            "GA1:B2M",
            "GA1L:B2M",
        ]

        for pair in combinations(chains, 2):

            # Get the chain objects
            c1, c2 = pair
            # What type of MHC are we pairing?
            the_type = ":".join(sorted([c1.chain_type, c2.chain_type], reverse=True))

            if the_type in acceptable_types:
                # Sort by chain type;
                p1 = pair[0] if c1.chain_type > c2.chain_type else pair[1]
                p2 = pair[1] if c2.chain_type < c1.chain_type else pair[0]

                try:
                    a1, a2 = (
                        p1[points[p1.chain_type][0]]["CA"],
                        p1[points[p1.chain_type][1]]["CA"],
                    )
                    a3, a4 = (
                        p2[points[p2.chain_type][0]]["CA"],
                        p2[points[p2.chain_type][1]]["CA"],
                    )
                    dist_array = [a3 - a1, a4 - a2, a4 - a1, a3 - a2]
                    constants = all(
                        [
                            dist_array[i] <= MHC_CUTOFF[the_type][i]
                            for i in range(len(dist_array))
                        ]
                    )
                    if constants:
                        pairings.append(pair)

                except KeyError:
                    continue

        return pairings

    def _get_sugar_fragments(self, sugar):
        """
        Get connected hetatoms to form sugar molecules.
        """
        # Make a sugar dictionary
        sugar = dict(list(zip([s.id for s in sugar], sugar)))

        # Get the connect records for the bonded atoms
        #        1 -  6         Record name      "CONECT"
        #        7 - 11          Integer          serial          Atom serial number
        #        12 - 16         Integer          serial          Serial number of bonded atom
        #        17 - 21         Integer          serial          Serial number of bonded atom
        #        22 - 26         Integer          serial          Serial number of bonded atom
        #        27 - 31         Integer          serial          Serial number of bonded atom
        connect_records = {}
        for c in [
            line.strip() for line in self.current_parser.trailer if "CONECT" in line
        ]:
            try:
                connect_records[int(c[6:11])] = []
            except IndexError:
                continue
            for b, e in [(11, 16), (16, 21), (21, 26), (26, 31)]:
                try:
                    if c[b:e].strip():
                        connect_records[int(c[6:11])].append(int(c[b:e]))
                    else:
                        break
                except IndexError:
                    break
                except ValueError:
                    self.warnings.write(
                        "Warning: unexpected CONECT record format %s" % c.strip()
                    )

        monomer_atoms = []
        polymers = []
        if connect_records:
            # Get the serial_numbers to residue id.
            atomid_to_resid = {}
            for r in sugar:
                for atom in sugar[r]:
                    atomid_to_resid[atom.serial_number] = sugar[r].id

            # Get the residue connections
            r_connections = {}
            for a in connect_records:
                if a in atomid_to_resid:
                    try:
                        r_connections[atomid_to_resid[a]].update(
                            [
                                atomid_to_resid[ai]
                                for ai in connect_records[a]
                                if ai in atomid_to_resid
                            ]
                        )
                    except KeyError:
                        r_connections[atomid_to_resid[a]] = set(
                            [
                                atomid_to_resid[ai]
                                for ai in connect_records[a]
                                if ai in atomid_to_resid
                            ]
                        )

            connected_sets = []
            for r in sorted(r_connections, key=lambda x: x[1]):
                added = 0
                for i in range(len(connected_sets)):
                    if connected_sets[i] & r_connections[r]:
                        connected_sets[i].update(r_connections[r])
                        added = 1
                        break
                if not added:
                    connected_sets.append(r_connections[r])

            n = 0
            for mol in connected_sets:
                if len(mol) > 1:
                    polymers.append(Fragment("sugar%d" % n))
                    for r in sorted(mol, key=lambda x: x[1]):
                        polymers[n].add(sugar[r])
                    n += 1
                else:
                    monomer_atoms += [atom for atom in sugar[list(mol)[0]]]

        else:
            for s in sugar:
                monomer_atoms += [atom for atom in sugar[s]]

        return polymers, monomer_atoms

    def _find_chain_hetatoms(self, chain):
        """
        Function for TCR and MHC chains to filter out hetatom records; this is to clean up the _match_units code below
        """
        hetatoms, sugars = [], []
        for residue in chain.get_unnumbered():
            # Ignore waters and non-standard amino acids
            if residue.id[0] == "W" or is_aa(residue, standard=False):
                continue
            if is_carbohydrate(residue):
                sugars.append(residue)
            else:
                hetatoms.extend(list(residue.get_atoms()))

        return hetatoms, sugars

    def _match_units(self, model, trchains, mhchains, agchains, crystal_contacts=[]):
        """
        Match MHC+Peptide chains to TCR chains.
        model is the current model - extract the TCRs from it (paired chains have been removed)
        trchains contains those TCR chains that have been unable to be paired to form TCRs
        agchains contains non-TCR chains that are potential antigens.

        Goal: Match TCR <-> MHC + peptide antigen.
        """
        # Get all T-cell receptor-like objects (TCR, TCRchain), and MHC-like objects.
        tcell_receptors = [h for h in model if isinstance(h, TCR)] + trchains.child_list
        mhc_complexes = [h for h in model if isinstance(h, MHC)] + mhchains.child_list

        # Initialise 5 dictionaries which carries a list of atoms per chain ID.
        antigen_atoms, cdr_atoms, mh_atoms, antigen_hetatoms, antigen_sugars = (
            defaultdict(list),
            defaultdict(list),
            defaultdict(list),
            defaultdict(list),
            defaultdict(list),
        )

        # Look through TCR and MHC and see if there are any weird hetatoms and sugars in the structure.
        for tr in tcell_receptors:
            for cdr in tr.get_CDRs():
                # Only get CDR3?
                if "3" not in cdr.id:
                    continue
                # only look at CA or CB atoms of the CDR; this is used later.
                cdr_atoms[tr.id] += [
                    atom
                    for atom in cdr.get_atoms()
                    if atom.id == "CB" or atom.id == "CA"
                ]

            # get TCR type and get chain's hetatoms accordingly
            if isinstance(tr, TCR) and tr.get_TCR_type() == "abTCR":
                beta_chain = tr.get_VB()
                alpha_chain = tr.get_VA()

                antigen_hetatoms[tr.VB], antigen_sugars[tr.VB] = (
                    self._find_chain_hetatoms(beta_chain)
                )
                antigen_hetatoms[tr.VA], antigen_sugars[tr.VA] = (
                    self._find_chain_hetatoms(alpha_chain)
                )

            elif isinstance(tr, TCR) and tr.get_TCR_type() == "gdTCR":
                delta_chain = tr.get_VD()
                gamma_chain = tr.get_VG()
                antigen_hetatoms[tr.VD], antigen_sugars[tr.VD] = (
                    self._find_chain_hetatoms(delta_chain)
                )
                antigen_hetatoms[tr.VG], antigen_sugars[tr.VG] = (
                    self._find_chain_hetatoms(gamma_chain)
                )

            elif isinstance(tr, TCR) and tr.get_TCR_type() == "dbTCR":
                beta_chain = tr.get_VB()
                delta_chain = tr.get_VD()
                antigen_hetatoms[tr.VB], antigen_sugars[tr.VB] = (
                    self._find_chain_hetatoms(beta_chain)
                )
                antigen_hetatoms[tr.VD], antigen_sugars[tr.VD] = (
                    self._find_chain_hetatoms(delta_chain)
                )

            # Unpaired TCR chain
            elif isinstance(tr, TCRchain):
                antigen_hetatoms[tr.id], antigen_sugars[tr.id] = (
                    self._find_chain_hetatoms(tr)
                )

        # Do the same for MHC.
        for mh in mhc_complexes:
            # Keep G domain atoms; Get the Helix region of MHC
            mh_atoms[mh.id] = [
                atom
                for atom in mh.get_atoms()
                if (atom.id == "CB" or atom.id == "CA") and atom.region == "Helix"
            ]
            if isinstance(mh, MHC) and mh.MHC_type == "MH1":
                MH1, B2M = mh.get_MH1(), mh.get_B2M()
                if MH1 is not None:
                    antigen_hetatoms[mh.MH1], antigen_sugars[mh.MH1] = (
                        self._find_chain_hetatoms(MH1)
                    )
                else:
                    GA1 = mh.get_GA1()
                    antigen_hetatoms[mh.GA1], antigen_sugars[mh.GA1] = (
                        self._find_chain_hetatoms(GA1)
                    )
                if B2M is not None:         # handle single chain MH1 case
                    antigen_hetatoms[mh.B2M], antigen_sugars[mh.B2M] = (
                        self._find_chain_hetatoms(B2M)
                    )

            elif isinstance(mh, MHC) and mh.MHC_type == "CD1":
                CD1, B2M = mh.get_CD1(), mh.get_B2M()
                if CD1 is not None:
                    antigen_hetatoms[mh.CD1], antigen_sugars[mh.CD1] = (
                        self._find_chain_hetatoms(CD1)
                    )
                if B2M is not None:
                    antigen_hetatoms[mh.B2M], antigen_sugars[mh.B2M] = (
                        self._find_chain_hetatoms(B2M)
                    )

            elif isinstance(mh, MHC) and mh.MHC_type == "MR1":
                MR1, B2M = mh.get_MR1(), mh.get_B2M()
                antigen_hetatoms[mh.MR1], antigen_sugars[mh.MR1] = (
                    self._find_chain_hetatoms(MR1)
                )
                antigen_hetatoms[mh.B2M], antigen_sugars[mh.B2M] = (
                    self._find_chain_hetatoms(B2M)
                )

            elif isinstance(mh, MHC) and mh.MHC_type == "MH2":
                GA, GB = mh.get_GA(), mh.get_GB()
                antigen_hetatoms[mh.GA], antigen_sugars[mh.GA] = (
                    self._find_chain_hetatoms(GA)
                )
                antigen_hetatoms[mh.GB], antigen_sugars[mh.GB] = (
                    self._find_chain_hetatoms(GB)
                )

            # Unpaired MHC chains -- if any, go here.
            elif isinstance(mh, MHCchain):
                antigen_hetatoms[mh.id], antigen_sugars[mh.id] = (
                    self._find_chain_hetatoms(mh)
                )

        for antigen in agchains:
            antigen_atoms[antigen.id] = [
                a
                for a in antigen.get_atoms()
                if a.parent.id[0] == " " or is_aa(a.parent)
            ]  # test ATOM records or amino acid HETATM records
            antigen_hetatoms[antigen.id] = [
                a
                for a in antigen.get_atoms()
                if a.parent.id[0].startswith("H") and not is_aa(a.parent)
            ]  # hetatm and not an amino acid

        # Problem here with carbohydrate molecules as units not recognised as polymers.
        # Have to use connect records to join them
        # Then consider them in the same way.
        sugars = []
        for chain_id in antigen_sugars:
            if antigen_sugars[chain_id]:
                polymers, monomer_atoms = self._get_sugar_fragments(
                    antigen_sugars[chain_id]
                )
                sugars += polymers
                antigen_hetatoms[chain_id] += monomer_atoms

        # We look through hetatms -- sometimes, hetatms can be associated to TR or MH chains, so we separate and whisk these out.
        # Protein/peptide entities override small molecules that are more likely to be buffer or cofactor molecules.
        # Get non-empty antigen hetatoms first
        non_empty_ag = [k for k in antigen_hetatoms if antigen_hetatoms[k]]

        # Pair proteins/peptides with TCR then MHC
        self._protein_peptide_pass(
            model, tcell_receptors, cdr_atoms, antigen_atoms, crystal_contacts
        )
        self._het_sugar_pass(
            tcell_receptors,
            cdr_atoms,
            non_empty_ag,
            antigen_hetatoms,
            sugars,
            distance=8.0,
        )

        # If a TCR does not have a detected MHC chain, then skip the remaining MHC-specific parsing bits.
        if not mhc_complexes:
            return

        # Have a very tight cutoff for MHCs that present het atoms (e.g. CD1 types)
        self._het_sugar_pass(
            mhc_complexes,
            mh_atoms,
            non_empty_ag,
            antigen_hetatoms,
            sugars,
            distance=3.5,
        )

        if antigen_atoms:
            self._protein_peptide_pass(
                model, mhc_complexes, mh_atoms, antigen_atoms, crystal_contacts
            )

        # Pair a TCR with an MHC and vice-versa; go through all possible combinations of TCR/MHC
        # We see if a CB/CA atom of the helix region of an MHC is within 8A of a TCR CDR loop's CB/CA atoms.
        # This is similar to the _protein_peptide_pass algorithm; we find the number of contacts between MHC and TCR,
        # and use the MHC with highest no. of contacts
        contact_freq = defaultdict(int)
        tr_mh_pairs = list(product(tcell_receptors, mhc_complexes))

        for tr, mh in tr_mh_pairs:
            ns = NeighborSearch(cdr_atoms[tr.id])
            for atom in mh_atoms[mh.id]:
                # This is a generous cutoff to be used for now.
                contacts = ns.search(atom.get_coord(), 8.0, level="R")
                for c in contacts:
                    contact_freq[(tr.id, mh.id)] += 1

        # Sort TR-MH pairs by number of contacts and then get the highest-frequency pairs
        sorted_contacts = sorted(
            list(contact_freq.items()), key=lambda z: z[1], reverse=True
        )
        paired_tr_mh = set()
        for pair, contacts in sorted_contacts:
            tr, mh = pair
            # If the TCR has already been paired, or if we know that the TCR and MHC are forming crystal contacts, move on.
            if tr in paired_tr_mh or (tr, mh) in crystal_contacts:
                continue
            model[tr]._add_mhc(model[mh])
            model[mh]._add_tcr(model[tr])
            paired_tr_mh.add(tr)

    def _protein_peptide_pass(
        self, model, complexes, receptor_atoms, antigen_atoms, crystal_contacts=[]
    ):
        """
        This is a generic method to process which proteins/peptides belong to a TCR or MHC. Needs testing.

        Args:
            complexes:       list of TCR/TCRchain objects or MHC/MHCchain objects
            receptor_atoms:  list of atom subset that will likely contact the antigen (e.g. cdr_atoms)
            antigen_atoms:   list of atoms in the antigen.
        """
        ns = NeighborSearch(
            [atom for chain in receptor_atoms for atom in receptor_atoms[chain]]
            + [atom for chain in antigen_atoms for atom in antigen_atoms[chain]]
        )
        contacts = [con for con in ns.search_all(8.0, "R")]
        contact_freq = defaultdict(lambda: defaultdict(int))

        # all_cpx_chains is a dictionary that has a TCR/MHC chain as a key and the ID of the TCR/MHC as value
        all_cpx_chains = dict()
        for cpx in complexes:
            cpx_ch = list(cpx.id)
            for c in cpx_ch:
                all_cpx_chains[c] = cpx.id

        # trids stores all paired/unpaired TR chains
        cpxids = set(all_cpx_chains.values())
        ags = set()

        for c in contacts:
            p1 = str(c[0].parent.id)  # get the chain id
            p2 = str(c[1].parent.id)

            # Reject cases where contacts are from the same chain, or the combination of chains is a TCR
            potential_contact = p1 + p2
            potential_contact2 = p2 + p1

            if (
                p1 == p2
                or potential_contact in contact_freq
                or potential_contact2 in contact_freq
            ):
                continue

            # If the potential contacting set of chains (p1+p2) is not a TR and p1 is a TR chain but p2 is NOT a TR chain, then p2 is an AG
            if (
                (potential_contact not in cpxids)
                and (p1 in all_cpx_chains)
                and (p2 not in all_cpx_chains)
            ):
                T = all_cpx_chains[p1]
                ag = p2
            # If the second set of potential contacting set of chains (p2+p1) is not a TR and p2 is a TR chain but p1 is NOT a TR chain,
            # then p1 is an AG
            elif (
                (potential_contact2 not in cpxids)
                and (p2 in all_cpx_chains)
                and (p1 not in all_cpx_chains)
            ):
                T = all_cpx_chains[p2]
                ag = p1
            else:
                continue

            # T is either the paired TCR id or an id of a single TCR chain
            contact_freq[T][ag] += 1
            ags.add(ag)

        # Iterate over the TR identifiers
        for cpx_id in cpxids:

            # If there are detected antigen contacts
            if contact_freq[cpx_id]:
                # Get the antigen
                ag = max(contact_freq[cpx_id], key=lambda x: contact_freq[cpx_id][x])

                if (cpx_id, ag) not in crystal_contacts:
                    model[cpx_id].antigen = (
                        []
                    )  # disregard smaller antigens if peptide or protein present.
                    model[cpx_id]._add_antigen(
                        model[ag]
                    )  # pair up an antigen to the TCR.

                    # Remove the antigen now, as it is paired up with a TR.
                    if ag in ags:
                        ags.remove(ag)

        # iterate over the remaining antigens to see if they are also bound.
        for ag in ags:
            cmax = 0
            for C in contact_freq:
                if ag in contact_freq[C] and (C, ag) not in crystal_contacts:
                    if contact_freq[C][ag] > cmax:
                        paired_cpx = C
                        cmax = contact_freq[C][ag]
            if cmax:
                if len(contact_freq) > 1:
                    self.warnings.write(
                        "Crystal Contact Warning: antigen %s has been paired with TCR %s"
                        % (str(ag), str(paired_cpx))
                    )
                    model[paired_cpx]._add_antigen(model[ag])
                else:
                    model[paired_cpx]._add_antigen(model[ag])

    def _het_sugar_pass(
        self,
        receptors,
        receptor_atoms,
        non_empty_ag,
        antigen_hetatoms,
        sugars,
        distance=8.0,
    ):
        """ """
        # Iterate through every possible pair of TR and hetatom chain
        for rec, antigen_het in product(receptors, non_empty_ag):
            # Initialise a NeighborSearch based on the atoms for a particualr chain of hetatoms
            ns = NeighborSearch(antigen_hetatoms[antigen_het])

            # Look through CDR atoms (CA/CB)
            for atom in receptor_atoms[rec.id]:
                # use 8.0A from the CDR CA/CB to the antigen. Using level = "R" returns Residue objects.
                contacts = ns.search(atom.get_coord(), distance, level="R")
                if contacts:
                    for contact in contacts:
                        # we assume that each contact residue is a single molecule (need to test its not just a residue)
                        if self._check_het_antigen(contact):
                            residue_type = get_res_type(contact)

                            if residue_type == "Hapten":
                                self.warnings.write(
                                    """Warning: Multiple hapten-antigen like molecules found in binding site -
                                    this needs attention as could be solvent/cofactor."""
                                )
                            if residue_type == "non-polymer":
                                contact.type = "Hapten"  # add a antigen type attribute to the residue
                                contact.get_type = (
                                    lambda: "Hapten"
                                )  # add a get antigen type method to the residue
                            elif residue_type == "nucleic-acid":
                                contact.type = "nucleic-acid"  # add a antigen type attribute to the residue
                                contact.get_type = (
                                    lambda: "nucleic-acid"
                                )  # add a get antigen type method to the residue
                            elif residue_type == "saccharide":
                                contact.type = "carbohydrate"  # add a antigen type attribute to the residue
                                contact.get_type = (
                                    lambda: "carbohydrate"
                                )  # add a get antigen type method to the residue
                            rec._add_antigen(contact)

        # Iterate through sugar fragments
        # for rec, sugar_fragment in product( receptors, sugars ):
        #     ns = NeighborSearch([atom for atom in sugar_fragment.get_atoms()])
        #     for atom in receptor_atoms[rec.id]:
        #         contacts = ns.search(atom.get_coord(), distance, level="R")
        #         if contacts:
        #             sugar_fragment.type = "carbohydrate" # add a antigen type attribute to the fragment
        #             sugar_fragment.get_type = lambda: "carbohydrate" # add a get antigen type method to the fragment
        #             rec._add_antigen(sugar_fragment)

    def _check_het_antigen(self, residue):
        """
        Method to perform checks on a potential hetatm residue.

        1. Check that it is not an amino acid - we don't want a modified residue to be found as a hapten.
        2. Check that the residue name is not a common buffer using high frequency residue codes.

        If we throw it out due to check 3 it will be reported to user.
        """

        # check 1
        # no amino acids
        if is_aa(residue, standard=False):
            return False

        # check 2
        # check for common buffers/unlikely haptens
        if is_common_buffer(residue):
            if not self.QUIET:
                self.warnings.write(
                    "Common molecule %s found in the binding site - not considered an antigen"
                    % residue.get_resname()
                )
            return False

        # add more checks as problems arise
        return True

    def _prenumbered(self, chain, prenumbering, ali_dict={}, n=0):
        """
        Method to deal with numbering supplied by the user. (or from the database)
        """

        if ali_dict:
            ali_dict = ali_dict[chain.id][n]

        annotation, chain_type = prenumbering[chain.id][n]

        try:
            sequence_list, sequence_str, warnings = extract_sequence(
                chain, return_warnings=True
            )
            numbering = align_numbering(annotation, sequence_list, ali_dict)
        except (
            AssertionError
        ):  # If the user has an alignment file generated before hetatoms included
            sequence_list, sequence_str, warnings = extract_sequence(
                chain, return_warnings=True, ignore_hets=True
            )
            numbering = align_numbering(annotation, sequence_list, ali_dict)
        self.warnings.log += warnings

        return numbering, chain_type


# class error_stream:
#     def __init__(self):
#         self.log = []

#     def __str__(self):
#         return "\n".join(self.log)

#     def __repr__(self):
#         return self.__str__()

#     def write(self, s):
#         self.log.append(str(s).strip("\n"))
