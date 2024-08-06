"""
Created on 9 May 2017
@author: leem
Modified version of the ABDB.AbPDB.Fragment class
"""

from .Entity import Entity
from Bio.PDB.PDBExceptions import PDBConstructionException


class Fragment(Entity):
    """
    A modified Entity class that can be thought of as a way of grouping children:
        e.g. TCR (TCR object) -> TCRchain (TCRchain object) -> Fragment CDRB3 (Fragment object)
                                                            -> Residue  B110  (Residue object)
    Does not modify the parent/child attributes of its children.
    For instance, one might define a fragment and add residues to it in order to visualise them.
    """

    def __init__(self, id):
        self._id = id
        Entity.__init__(self, id)
        self.level = "F"

    def __repr__(self):
        if hasattr(self, "chain_type"):
            return "<Fragment %s TCRchain: %s>" % (self.id, self.parent.parent.id)
        else:
            return "<Fragment %s>" % self.id

    def add(self, entity):
        "Add a child to the Entity."
        entity_id = entity.get_id()
        if self.has_id(entity_id):
            raise PDBConstructionException("%s defined twice" % str(entity_id))

        # parent of child is not changed
        self.child_list.append(entity)
        self.child_dict[entity_id] = entity

    def insert(self, pos, entity):
        "Add a child to the Entity at a specified position."
        entity_id = entity.get_id()
        if self.has_id(entity_id):
            raise PDBConstructionException("%s defined twice" % str(entity_id))

        # parent of child is not changed
        self.child_list[pos:pos] = [entity]
        self.child_dict[entity_id] = entity

    def get_residues(self):
        for residue in self:
            yield residue

    def get_atoms(self):
        for residue in self.get_residues():
            for atom in residue:
                yield atom
