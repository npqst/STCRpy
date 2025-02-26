"""
Created on 30 Apr 2016
@author: leem

Based on the ABchain class from @dunbar

"""

from .utils.region_definitions import get_region
from Bio.PDB.Chain import Chain
from .Entity import Entity
from .Fragment import Fragment
from Bio import SeqUtils


class MHCchain(Chain, Entity):
    """
    A class to hold an MHC chain
    """

    def __init__(self, id):
        Chain.__init__(self, id)
        Entity.__init__(self, id)
        self.level = "C"
        self.antigen = []
        self.unnumbered = []
        self.sequence = {}
        self.residue_order = {}
        self.engineered = False
        self.tcr = []

    def __repr__(self):
        return "<MHCchain %s type: %s>" % (self.id, self.chain_type)

    def _add_antigen(self, antigen=None):
        if antigen not in self.antigen:
            self.antigen.append(antigen)

    def _add_tcr(self, tcr):
        self.tcr.append(tcr)

    def get_TCR(self):
        return self.tcr

    def analyse(self, chain_type):
        self.set_chain_type(chain_type)
        # self._init_fragments()
        self.annotate_children()
        self.set_sequence()

    def set_chain_type(self, chain_type):
        """
        Set the MHC's chain type
        """
        self.chain_type = chain_type

    def set_engineered(self, engineered):
        if engineered:
            self.engineered = True
        else:
            self.engineered = False

    def set_sequence(self):
        i = 0
        for residue in self:
            if (
                residue.get_resname().capitalize()
                in SeqUtils.IUPACData.protein_letters_3to1
            ):
                resname = SeqUtils.IUPACData.protein_letters_3to1[
                    residue.get_resname().capitalize()
                ]  # change this to use our chemical components.
            else:
                # skip the residue if the code is not recognised - e.g. UNK
                continue
            hetflag, resseq, icode = residue.get_id()
            self.sequence[(self.chain_type + str(resseq) + str(icode)).strip()] = (
                resname
            )
            self.residue_order[(self.chain_type + str(resseq) + str(icode)).strip()] = i
            i += 1

    def add_unnumbered(self, residue):
        self.unnumbered.append(residue.id)

    def _get_region(self, residue):
        region = ""
        if hasattr(residue, "imgt_numbered") and residue.imgt_numbered:
            # Call the get_region function from ..Annotate.region_definitions module
            region = get_region((residue.id[1], residue.id[2]), residue.chain_type)
            return region
        return "?"

    def annotate_children(self):
        for residue in self:
            residue.chain_type = self.chain_type
            residue.region = self._get_region(residue)
            for atom in residue:
                atom.chain_type = self.chain_type
                atom.region = residue.region
            # if residue.region != "?":
            #    self.fragments.child_dict[residue.region].add(residue)

    def _init_fragments(self):
        self.fragments = Entity("Fragments")
        self.fragments.set_parent(self)
        for region in regions[self.chain_type]:
            self.fragments.add(Fragment(region))

    def is_bound(self):
        """
        Check whether there is an antigen bound to the antibody chain
        """
        if self.get_antigen():
            return True
        else:
            return False

    def get_antigen(self):
        return self.antigen

    def get_fragments(self):
        for f in self.fragments:
            yield f

    def get_sequence(self, type=dict):
        if not self.sequence:
            self.set_sequence()
        if type is dict:
            return self.sequence
        else:
            ordered = sorted(
                list(self.sequence.items()), key=lambda x: self.residue_order[x[0]]
            )
            if type is str:
                return "".join([r[1] for r in ordered])
            else:
                return ordered

    def get_unnumbered(self):
        for r in self.unnumbered:
            yield self.child_dict[r]

    def get_chains(self):  # implemented to retain interface with MHC.get_chains()
        for c in [self]:
            yield c

    def get_allele_assignments(self):
        return self.xtra["genetic_origin"]
