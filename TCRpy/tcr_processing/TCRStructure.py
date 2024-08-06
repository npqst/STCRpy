"""
Created on 10 May 2017
@author: leem
Based on the ABDB.AbPDB.AntibodyStructure class.
"""

from Bio import SeqUtils
from .Entity import Entity
from .TCR import TCR
from .MHC import MHC


class TCRStructure(Entity):
    """
    The TCRStructure class contains a collection of models
    """

    def __init__(self, identifier):
        self.level = "TS"
        Entity.__init__(self, identifier)
        self.header = {}

    def __repr__(self):
        return "<Structure id=%s>" % self.get_id()

    def _sort(self, m1, m2):
        """Sort models.

        This sorting function sorts the Model instances in the Structure instance.
        The sorting is done based on the model id, which is a simple int that
        reflects the order of the models in the PDB file.

        Arguments:
        o m1, m2 - Model instances
        """
        return (m1.get_id() > m2.get_id()) - (m1.get_id() < m2.get_id())

    def _set_numbering_scheme(self, scheme=None):
        """
        Set the numbering scheme used.
        """
        self.numbering_scheme = scheme

    # Public
    def set_header(self, header):
        """
        Set the header as the parsed header dictionary from biopython
        """
        self.header = header

    def get_header(self):
        return self.header

    def get_models(self):
        for m in self:
            yield m

    def get_holders(self):
        for m in self.get_models():
            for h in m:
                yield h

    def get_TCRs(self):
        """
        Get any instance of the TCR object.
        Hierarchy:
            TCRStructure
               |
               |______ TCR
               |
               |______ MHC
        """
        for h in self.get_holders():
            if isinstance(h, TCR):
                yield h

    def get_TCRchains(self):
        """Gets all TCR chains"""
        for h in self.get_holders():
            if h.id == "TCRchain":
                for c in h:
                    yield c
            elif isinstance(h, TCR):
                for c in h:
                    yield c

    def get_MHCs(self):
        """
        Get any instance of the MHC object.
        Hierarchy:
            TCRStructure
               |
               |______ TCR
               |
               |______ MHC
        """
        for h in self.get_holders():
            if isinstance(h, MHC):
                yield h

    def get_antigens(self):
        """
        This gets the 'antigen' chains in the structure,
        that have been assigned to a TCR or an MHC.
        """
        antigens = set([])
        for h in self.get_holders():
            if isinstance(h, MHC) or isinstance(h, TCR) or h.id == "TCRchain":
                for c in h.antigen:
                    if c not in antigens:
                        antigens = antigens.union(set([c]))
                        yield c

    def get_unpaired_TCRchains(self):
        """
        This gets the TCR chains that are not paired
        """
        for h in self.get_holders():
            if h.id == "TCRchain":
                for c in h:
                    yield c

    def get_chains(self):
        for h in self.get_holders():
            for c in h:
                yield c

    def get_residues(self):
        for c in self.get_chains():
            for r in c:
                yield r

    def get_atoms(self):
        for r in self.get_residues():
            for a in r:
                yield a

    def get_seq(self, model=0):
        seq = ""
        for c in self[model]:
            for r in c.get_residues():
                # Skip over water molecules
                if r.resname == "HOH":
                    continue
                seq += SeqUtils.IUPACData.protein_letters_3to1[r.resname]
            seq += "/"

        return seq[:-1]
