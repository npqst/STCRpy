"""
Created on 30 Apr 2016

@author: leem, based on work by dunbar

The MHC class. This is similar to the Fab class.

"""

from .Entity import Entity


class MHC(Entity):
    """
    MHC class.
    Holds paired MHC domains.
    """

    def __init__(self, c1, c2):
        if hasattr(c1, "chain_type"):
            Entity.__init__(self, c1.id + c2.id)
        else:
            Entity.__init__(self, c2.id + c1.id)

        self.level = "H"
        self._add_domain(c1)
        self._add_domain(c2)
        self._set_MHC_type()
        self.child_list = sorted(self.child_list, key=lambda x: x.id)
        self.antigen = []
        self.tcr = []
        self.engineered = False

    def _add_antigen(self, antigen=None):
        if antigen not in self.antigen:
            self.antigen.append(antigen)

    def _add_tcr(self, tcr=None):
        self.tcr.append(tcr)

    def get_TCR(self):
        return self.tcr

    def get_antigen(self):
        """
        Return a list of bound antigens.
        If the antigen has more than one chain, those in contact with the antibody will be returned.
        """
        return self.antigen

    def is_bound(self):
        """
        Check whether there is an antigen bound to the antibody fab
        """
        if self.get_antigen():
            return True
        else:
            return False

    def get_chains(self):
        for c in self:
            yield c

    def get_residues(self):
        for c in self.get_chains():
            for r in c:
                yield r

    def get_atoms(self):
        for r in self.get_residues():
            for a in r:
                yield a

    def get_MHC_type(self):
        if hasattr(self, "MHC_type"):
            return self.MHC_type

    def get_allele_assignments(self):
        return {c.id: c.get_allele_assignments() for c in self.get_chains()}


class MH1(MHC):
    """
    Type 1 MHC class.
    Holds paired MHC domains.
    """

    def __repr__(self):
        if self.MHC_type == "MH1":
            return "<%s %s%s GA1/GA2=%s; B2M=%s>" % (
                self.MHC_type,
                self.MH1,
                self.B2M,
                self.MH1,
                self.B2M,
            )
        else:
            return "<GA1/GA2 %s%s GA1=%s; GA2=%s>" % (
                self.GA1,
                self.GA2,
                self.GA1,
                self.GA2,
            )

    def _set_MHC_type(self):
        if hasattr(self, "MH1"):
            self.MHC_type = "MH1"
        elif hasattr(self, "GA1") or hasattr(self, "GA2"):
            self.MHC_type = "MH1"
        else:
            self.MHC_type = ""

    def _add_domain(self, chain):
        if chain.chain_type == "MH1":
            self.MH1 = chain.id
            self.GA1 = chain.id
            self.GA2 = chain.id
        elif chain.chain_type == "GA1":
            self.MH1 = chain.id
            self.GA1 = chain.id
        elif chain.chain_type == "GA2":
            self.MH1 = chain.id
            self.GA2 = chain.id
        elif chain.chain_type == "B2M":
            self.B2M = chain.id

        # Add the chain as a child of this entity.
        self.add(chain)

    def get_alpha(self):
        for MH1_domain in set(["MH1", "GA1", "GA2"]):
            if hasattr(self, MH1_domain):
                return self.child_dict[getattr(self, MH1_domain)]

    def get_MH1(self):
        if hasattr(self, "MH1"):
            return self.child_dict[self.MH1]

    def get_GA1(self):
        if hasattr(self, "GA1"):
            return self.child_dict[self.GA1]
        else:
            return self.get_MH1()

    def get_GA2(self):
        if hasattr(self, "GA2"):
            return self.child_dict[self.GA2]
        else:
            return self.get_MH1()

    def get_B2M(self):
        if hasattr(self, "B2M"):
            return self.child_dict[self.B2M]


class MH2(MHC):
    """
    Type 2 MHC class.
    Holds paired MHC domains.
    """

    def __repr__(self):
        if self.MHC_type == "MH2":
            return "<%s %s%s GA=%s; GB=%s>" % (
                self.MHC_type,
                self.GA,
                self.GB,
                self.GA,
                self.GB,
            )
        else:
            return "<GA/GB %s%s GA=%s; GB=%s>" % (self.GA, self.GB, self.GA, self.GB)

    def _set_MHC_type(self):
        if hasattr(self, "GA") and hasattr(self, "GB"):
            self.MHC_type = "MH2"
        else:
            self.MHC_type = ""

    def _add_domain(self, chain):
        if chain.chain_type == "GA":
            self.GA = chain.id
        elif chain.chain_type == "GB":
            self.GB = chain.id

        # Add the chain as a child of this entity.
        self.add(chain)

    def get_GA(self):
        if hasattr(self, "GA"):
            return self.child_dict[self.GA]

    def get_GB(self):
        if hasattr(self, "GB"):
            return self.child_dict[self.GB]


class CD1(MHC):
    """
    CD1 class.
    Holds paired CD1/B2M domains.
    """

    def __repr__(self):
        if self.MHC_type == "CD1":
            return "<%s %s%s GA1L/GA2L=%s; B2M=%s>" % (
                self.MHC_type,
                self.CD1,
                self.B2M,
                self.CD1,
                self.B2M,
            )
        else:
            return "<GA1L/GA2L %s%s GA1L=%s; GA2L=%s>" % (
                self.GA1L,
                self.GA2L,
                self.GA1L,
                self.GA2L,
            )

    def _set_MHC_type(self):
        if hasattr(self, "CD1"):
            self.MHC_type = "CD1"
        elif hasattr(self, "GA1L") and hasattr(self, "GA2L"):
            self.MHC_type = "CD1"
        elif hasattr(self, "GA1L") and hasattr(self, "B2M"):
            self.MHC_type = "CD1"
        else:
            self.MHC_type = ""

    def _add_domain(self, chain):
        if chain.chain_type == "CD1":
            self.MHC_type = "CD1"
            self.CD1 = chain.id
            self.GA1L = chain.id
            self.GA2L = chain.id
        elif chain.chain_type == "GA1L":
            self.CD1 = chain.id
            self.GA1L = chain.id
        elif chain.chain_type == "GA2L":
            self.GA2L = chain.id
        elif chain.chain_type == "B2M":
            self.B2M = chain.id

        # Add the chain as a child of this entity.
        self.add(chain)

    def get_CD1(self):
        if hasattr(self, "CD1"):
            return self.child_dict[self.CD1]

    def get_B2M(self):
        if hasattr(self, "B2M"):
            return self.child_dict[self.B2M]


class MR1(MHC):
    """
    MR1 class.
    Holds paired MR1/B2M domains.
    """

    def __repr__(self):
        if self.MHC_type == "MR1":
            return "<%s %s%s GA1L/GA2L=%s; B2M=%s>" % (
                self.MHC_type,
                self.MR1,
                self.B2M,
                self.MR1,
                self.B2M,
            )
        else:
            return "<GA1L/GA2L %s%s GA1L=%s; GA2L=%s>" % (
                self.GA1L,
                self.GA2L,
                self.GA1L,
                self.GA2L,
            )

    def _set_MHC_type(self):
        if hasattr(self, "MR1"):
            self.MHC_type = "MR1"
        elif hasattr(self, "GA1L") and hasattr(self, "GA2L"):
            self.MHC_type = "MR1"
        else:
            self.MHC_type = ""

    def _add_domain(self, chain):
        if chain.chain_type == "MR1":
            self.MHC_type = "MR1"
            self.MR1 = chain.id
            self.GA1L = chain.id
            self.GA2L = chain.id
        elif chain.chain_type == "GA1L":
            self.GA1L = chain.id
        elif chain.chain_type == "GA2L":
            self.GA2L = chain.id
        elif chain.chain_type == "B2M":
            self.B2M = chain.id

        # Add the chain as a child of this entity.
        self.add(chain)

    def get_MR1(self):
        if hasattr(self, "MR1"):
            return self.child_dict[self.MR1]

    def get_B2M(self):
        if hasattr(self, "B2M"):
            return self.child_dict[self.B2M]


class scMH1(MHC):
    """
    Type 1 MHC class.
    Holds single chain MHC domains for Class I MHC if the identiifed chain
    is the double alpha helix, ie. MH1 without B2M, with exception for GA1.
    """

    def __init__(self, c1):
        assert c1.chain_type in [
            "GA1",
            "GA2",
            "MH1",
        ], f"Chain {c1} with can not form a single chain MHC class I."
        Entity.__init__(self, c1.id)

        self.level = "H"
        self._add_domain(c1)
        self._set_MHC_type()
        self.child_list = sorted(self.child_list, key=lambda x: x.id)
        self.antigen = []
        self.tcr = []
        self.engineered = False

    def __repr__(self):
        if self.MHC_type == "MH1":
            return "<%s %s GA1/GA2=%s>" % (
                self.MHC_type,
                self.MH1,
                self.MH1,
            )
        else:
            return "<GA1/GA2 %s GA1/GA2=%s>" % (
                self.GA1,
                self.GA1,
            )

    def _set_MHC_type(self):
        if hasattr(self, "MH1") or hasattr(self, "GA1") or hasattr(self, "GA2"):
            self.MHC_type = "MH1"
        else:
            self.MHC_type = ""

    def _add_domain(self, chain):
        if chain.chain_type in ["MH1", "GA1", "GA2"]:
            self.MH1 = chain.id
            self.GA1 = chain.id
            self.GA2 = chain.id

        # Add the chain as a child of this entity.
        self.add(chain)

    def get_alpha(self):
        for MH1_domain in set(["MH1", "GA1", "GA2"]):
            if hasattr(self, MH1_domain):
                return self.child_dict[getattr(self, MH1_domain)]

    def get_MH1(self):
        if hasattr(self, "MH1"):
            return self.child_dict[self.MH1]

    def get_GA1(self):
        if hasattr(self, "GA1"):
            return self.child_dict[self.GA1]
        else:
            return self.get_MH1()

    def get_GA2(self):
        if hasattr(self, "GA2"):
            return self.child_dict[self.GA2]
        else:
            return self.get_MH1()

    def get_B2M(self):
        return None


class scCD1(MHC):
    """
    Type 1 MHC class.
    Holds single chain MHC domains of type CD1 for Class I MHC if the identiifed chain
    is the double alpha helix, ie. CD1 without B2M.
    """

    def __init__(self, c1):
        assert c1.chain_type in [
            "GA1L",
            "CD1",
        ], f"Chain {c1} with can not form a single chain MHC class I."
        Entity.__init__(self, c1.id)

        self.level = "H"
        self._add_domain(c1)
        self._set_MHC_type()
        self.child_list = sorted(self.child_list, key=lambda x: x.id)
        self.antigen = []
        self.tcr = []
        self.engineered = False

    def __repr__(self):
        if self.MHC_type == "CD1":
            return "<%s %s GA1L=%s>" % (
                self.MHC_type,
                self.CD1,
                self.CD1,
            )
        else:
            return "<GA1L %s GA1L=%s>" % (
                self.GA1L,
                self.GA1L,
            )

    def _set_MHC_type(self):
        if hasattr(self, "CD1") or hasattr(self, "GA1L"):
            self.MHC_type = "CD1"
        else:
            self.MHC_type = ""

    def _add_domain(self, chain):
        if chain.chain_type in ["CD1", "GA1L"]:
            self.CD1 = chain.id
            self.GA1L = chain.id

        # Add the chain as a child of this entity.
        self.add(chain)

    def get_CD1(self):
        if hasattr(self, "CD1"):
            return self.child_dict[self.CD1]

    def get_GA1L(self):
        if hasattr(self, "GA1L"):
            return self.child_dict[self.GA1L]
        else:
            return self.get_CD1()

    def get_B2M(self):
        return None


class scMH2(MHC):
    """
    Single chain MHC class 2.
    Holds single GA or GB chain.
    Usually this will only occur if ANARCI has not been identified one of the two chains correctly.
    """

    def __init__(self, c1):
        assert c1.chain_type in [
            "GA",
            "GB",
        ], f"Chain {c1} with can not form a single chain MHC class I."
        Entity.__init__(self, c1.id)

        self.level = "H"
        self._add_domain(c1)
        self._set_MHC_type()
        self.child_list = sorted(self.child_list, key=lambda x: x.id)
        self.antigen = []
        self.tcr = []
        self.engineered = False

    def __repr__(self):
        if self.MHC_type == "MH2":
            if hasattr(self, "GA"):
                return "<%s %s GA=%s>" % (
                    self.MHC_type,
                    self.GA,
                    self.GA,
                )
            elif hasattr(self, "GB"):
                return "<%s %s GB=%s>" % (
                    self.MHC_type,
                    self.GB,
                    self.GB,
                )

        else:
            if hasattr(self, "GA"):
                return "<GA %s GA=%s>" % (self.GA, self.GA)
            elif hasattr(self, "GB"):
                return "<GB %s GB=%s>" % (self.GB, self.GB)

    def _set_MHC_type(self):
        if hasattr(self, "GA") and hasattr(self, "GB"):
            self.MHC_type = "MH2"
        else:
            self.MHC_type = ""

    def _add_domain(self, chain):
        if chain.chain_type == "GA":
            self.GA = chain.id
        elif chain.chain_type == "GB":
            self.GB = chain.id

        # Add the chain as a child of this entity.
        self.add(chain)

    def get_GA(self):
        if hasattr(self, "GA"):
            return self.child_dict[self.GA]

    def get_GB(self):
        if hasattr(self, "GB"):
            return self.child_dict[self.GB]
