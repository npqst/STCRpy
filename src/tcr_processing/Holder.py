"""
Created on 9 May 2017
@author: leem

A generic holder class that can be used to contain individual chains, etc.

"""

from .Entity import Entity


class Holder(Entity):
    def __init__(self, identifier):
        Entity.__init__(self, identifier)
        self.level = "H"

    def __repr__(self):
        if len(self.child_list):
            return "<Holder %s chains: %s>" % (
                self.id,
                ",".join([child.id for child in self]),
            )
        else:
            return "<Holder %s chains: None>" % (self.id)
