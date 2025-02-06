"""
Created on 9 May 2017
@author: leem

Based on the ABDB.AbPDB.Model class.
"""

import Bio


class Model(Bio.PDB.Model.Model, Bio.PDB.Entity.Entity):
    """
    Override to use our Entity

    @change: __getitem__ changed so that single chains can be called as well as holder object from a model.
            e.g. s[0]["B"] and s[0]["BA"] gets the B chain and the BA tcr respectively.
    """

    def __init__(self, identifier, serial_num=None):
        Bio.PDB.Model.Model.__init__(self, identifier, serial_num)
        Bio.PDB.Entity.Entity.__init__(self, identifier)

    def __getitem__(self, identifier):
        "Return the child with given identifier."
        if identifier in self.child_dict:
            return self.child_dict[identifier]
        else:
            # Allow a single chain to be called from a model.
            for child in self:
                if (
                    child.id in self.child_dict
                    and identifier in self.child_dict[child.id].child_dict
                ):
                    return self.child_dict[child.id].child_dict[identifier]
            raise KeyError(
                identifier
            )  # raise error if chain belonging to identifier is not found
