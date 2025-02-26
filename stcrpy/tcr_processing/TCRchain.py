from Bio.PDB import Chain
from Bio import SeqUtils
from .utils.region_definitions import get_region
from .Entity import Entity
from .Fragment import Fragment

regions = {
    "B": ["fwb1", "cdrb1", "fwb2", "cdrb2", "fwb3", "cdrb3", "fwb4"],
    "A": ["fwa1", "cdra1", "fwa2", "cdra2", "fwa3", "cdra3", "fwa4"],
    "D": ["fwd1", "cdrd1", "fwd2", "cdrd2", "fwd3", "cdrd3", "fwd4"],
    "G": ["fwg1", "cdrg1", "fwg2", "cdrg2", "fwg3", "cdrg3", "fwg4"],
}


class TCRchain(Chain.Chain, Entity):
    """
    A class to hold a TCR chain.
    """

    def __init__(self, identifier):
        Chain.Chain.__init__(self, identifier)
        Entity.__init__(self, identifier)
        self.level = "C"
        self.mhc = []
        self.antigen = []
        self.unnumbered = []
        self.sequence = {}
        self.residue_order = {}
        self.engineered = False

    def __repr__(self):
        return "<TCRchain %s type: %s>" % (self.id, self.chain_type)

    def _add_mhc(self, mhc=None):
        self.mhc.append(mhc)

    def _add_antigen(self, antigen=None):
        if antigen not in self.antigen:
            self.antigen.append(antigen)

    def is_bound(self):
        """
        Check whether there is an antigen bound to the TCR
        """
        if self.get_antigen():
            return True
        else:
            return False

    def analyse(self, chain_type):
        self.set_chain_type(chain_type)
        self._init_fragments()
        self.annotate_children()
        self.set_sequence()

    def set_chain_type(self, chain_type):
        """
        Set the chain type to B, A, D, or G
        """
        self.chain_type = chain_type

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

    def set_engineered(self, engineered):
        if engineered:
            self.engineered = True
        else:
            self.engineered = False

    def add_unnumbered(self, residue):
        self.unnumbered.append(residue.id)

    def _get_region(self, residue):
        region = ""
        if hasattr(residue, "imgt_numbered") and residue.imgt_numbered:
            region = get_region((residue.id[1], residue.id[2]), self.chain_type)
            return region
        return "?"

    def annotate_children(self):
        for residue in self:
            residue.chain_type = self.chain_type
            residue.region = self._get_region(residue)
            for atom in residue:
                atom.chain_type = self.chain_type
                atom.region = residue.region

            if residue.region != "?":
                self.fragments.child_dict[residue.region].add(residue)

    def _init_fragments(self):
        self.fragments = Entity("Fragments")
        self.fragments.set_parent(self)
        for region in regions[self.chain_type]:
            self.fragments.add(Fragment(region))

    def is_engineered(self):
        return self.engineered

    def get_MHC(self):
        return self.mhc

    def get_antigen(self):
        return self.antigen

    def get_fragments(self):
        for f in self.fragments:
            yield f

    def get_CDRs(self):
        for f in self.fragments:
            if f.id.lower()[:3] == "cdr":
                yield f

    def get_frameworks(self):
        """
        Obtain framework regions from a TCRChain object.
        """
        for f in self.get_fragments():
            if "fw" in f.id:
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

    def get_germline_assignments(self):
        return self.xtra["genetic_origin"]
