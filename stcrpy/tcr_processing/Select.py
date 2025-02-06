"""
Select.py
Created on 9 May 2017
@author: leem

These are selection classes for the save method of the TcrPDB entity
They are based on the ABDB.AbPDB.Select and Bio.PDB.PDBIO Selection classes

"""

from .utils.constants import BACKBONE_CB
from .utils.common import *


class select_all(object):
    """
    Default selection (everything) during writing - can be used as base class
    to implement selective output. This selects which entities will be written out.
    """

    def __repr__(self):
        return "<Select all>"

    def accept(self, ob):
        if ob.level == "A":
            return self.accept_atom(ob)
        elif ob.level == "R":
            return self.accept_residue(ob)
        elif ob.level == "C":
            return self.accept_chain(ob)
        elif ob.level == "F":
            return self.accept_fragment(ob)
        elif ob.level == "H":
            return self.accept_holder(ob)
        elif ob.level == "M":
            return self.accept_model(ob)

    def accept_model(self, model):
        """
        Overload this to reject models for output.
        """
        return 1

    def accept_holder(self, model):
        """
        Overload this to reject holders for output. (TCRs, TCRchains-holder, MHCchains-holder, AGchains-holder)
        """
        return 1

    def accept_chain(self, chain):
        """
        Overload this to reject chains for output.
        """
        return 1

    def accept_fragment(self, fragment):
        """
        Overload this to reject residues for output.
        """
        return 1

    def accept_residue(self, residue):
        """
        Overload this to reject residues for output.
        """
        return 1

    def accept_atom(self, atom):
        """
        Overload this to reject atoms for output.
        """
        return 1


class variable_only(select_all):
    """
    Select the variable region(s) of the structure.
    """

    def __repr__(self):
        return "<variable_only>"

    def accept_holder(self, holder):
        """
        Overload this to reject holders for output.
        """
        # Accept an abTCR or a gdTCR
        if (
            (hasattr(holder, "VB") and hasattr(holder, "VA"))
            or (hasattr(holder, "VD") and hasattr(holder, "VG"))
            or (hasattr(holder, "VB") and hasattr(holder, "VD"))
        ):
            return 1
        else:
            return 0

    def accept_residue(self, residue):
        """
        Overload this to reject residues for output.
        """
        if hasattr(residue, "region") and residue.region != "?":
            return 1
        else:
            return 0


class cdr3(variable_only):
    """
    Select only CDR3.
    """

    def __repr__(self):
        return "<CDR3>"

    def accept_residue(self, residue):
        if "cdr" in residue.region and "3" in residue.region:
            return 1
        else:
            return 0


class backbone(select_all):
    """
    Select only backbone (no side chains) atoms in the structure.
    Backbone defined as "C","CA","N","CB" and "O" atom identifiers in amino acid (pdb notation)
    """

    def __repr__(self):
        return "<backbone>"

    def accept_atom(self, atom):
        if atom.id in BACKBONE_CB:
            return 1
        else:
            return 0


class fv_only_backbone(variable_only, backbone):
    """
    Select the backbone atoms of the variable region.
    Example of combining selection classes.
    """

    def __repr__(self):
        return "<variable only backbone>"
