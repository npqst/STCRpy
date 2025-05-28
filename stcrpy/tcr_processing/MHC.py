"""
Created on 30 Apr 2016

@author: leem, based on work by dunbar

The MHC class. This is similar to the Fab class.

"""

import warnings

from Bio import BiopythonWarning

from .Entity import Entity
from .MHCchain import MHCchain
from .utils.region_definitions import IMGT_MH1_ABD, IMGT_MH2_ABD


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

    def copy(self, copy_siblings = True ):
        """
        Return a copy of the MHC object. This returns a shallow copy of the MHC object.
        If the copy_siblings flag is set to True, the antigen and TCR objects will also be copied. Warning - if the copy_siblings flag is set to False, the antigen and TCR objects will not be copied, and the reference will still point to the same TCR and antigen objects as the original.

        copy_siblings: Whether to copy sibling entities (ie. TCR and Antigen objects). Default True. 

        """
        shallow = super().copy()
        if copy_siblings:
            shallow.antigen = [a.copy(copy_siblings=False) for a in self.get_antigen()]
            shallow.tcr = [t.copy(copy_siblings=False) for t in self.get_TCR()]
            for t in shallow.tcr:
                t.MHC = [shallow if m.id == shallow.id else m.copy(copy_siblings=False) for m in t.get_MHC()]
                t.antigen = [ag.copy(copy_siblings=False) if ag.id not in [a.id for a in shallow.antigen] else [a for a in shallow.antigen if a.id==ag.id][0] for ag in t.antigen]
        return shallow

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

    def crop(self, *args, **kwargs):
        """Raises NotImplementedError."""
        raise NotImplementedError()

    def standardise_chain_names():
        """Raises NotImplementedError."""
        raise NotImplementedError()


class MH1(MHC):
    """
    Type 1 MHC class.
    Holds paired MHC domains.
    """

    def __repr__(self):
        if self.MHC_type == "MH1":
            return "<%s %s%s GA1/GA2=%s; B2M=%s>" % (
                self.MHC_type,
                self.MH1 if self.MH1 else '',
                self.B2M if self.B2M else '',
                self.MH1,
                self.B2M,
            )
        else:
            return "<GA1/GA2 %s%s GA1=%s; GA2=%s>" % (
                self.GA1 if self.GA1 else '',
                self.GA2 if self.GA1 else '',
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

    def crop(self, *, remove_het_atoms: bool = True) -> None:
        """Crop to antigen binding domain.

        This method mutates the MH1 object.

        Args:
            remove_het_atoms: remove het atoms from structure as well

        """
        alpha_chain = self.get_alpha()
        new_alpha_chain = MHCchain(alpha_chain.id)

        for residue in alpha_chain:
            if residue.id[1] in IMGT_MH1_ABD or (not remove_het_atoms and residue.id[0] != ' '):
                new_alpha_chain.add(residue.copy())

        new_alpha_chain.analyse(alpha_chain.chain_type)

        del self[alpha_chain.id]
        self.add(new_alpha_chain)

        if hasattr(self, 'B2M'):
            del self[self.B2M]
            self.B2M = None

    def standardise_chain_names(self) -> None:
        """Standardise MHC chain name to A and B2M chain name to B."""
        new_id = []
        new_child_dict = {}

        for MH1_domain in set(['MH1', 'GA1', 'GA2']):
            if hasattr(self, MH1_domain):
                new_child_dict['A'] = self.child_dict[getattr(self, MH1_domain)]
                setattr(self, MH1_domain, 'A')
                new_id.append('A')
                break

        if hasattr(self, 'B2M'):
            new_child_dict['B'] = self.child_dict[self.B2M]
            self.B2M = 'B'
            new_id.append('B')

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)

            for chain_id, chain in new_child_dict.items():
                chain.id = chain_id

        self.child_dict = new_child_dict

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            self.id = ''.join(new_id)


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

    def crop(self, *, remove_het_atoms: bool = True) -> None:
        """Crop to antigen binding domain.

        This method mutates the MH2 object.

        Args:
            remove_het_atoms: remove het atoms from structure as well

        """
        new_child_dict = {}

        for chain in self:
            new_chain = MHCchain(chain.id)

            for residue in chain:
                if residue.id[1] in IMGT_MH2_ABD or (not remove_het_atoms and residue.id[0] != ' '):
                    new_chain.add(residue.copy())

            new_chain.analyse(chain.chain_type)
            new_child_dict[new_chain.id] = new_chain

        for chain_id in new_child_dict:
            del self[chain_id]

        for new_chain in new_child_dict.values():
            self.add(new_chain)

    def standardise_chain_names(self) -> None:
        """Standardise MHC chain 1 name to A and MHC chain 2 name to B."""
        new_id = []
        new_child_dict = {}

        if hasattr(self, 'GA'):
            new_child_dict['A'] = self.child_dict[self.GA]
            self.GA = 'A'
            new_id.append('A')

        if hasattr(self, 'GB'):
            new_child_dict['B'] = self.child_dict[self.GB]
            self.GB = 'B'
            new_id.append('B')

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)

            for chain_id, chain in new_child_dict.items():
                chain.id = chain_id

        self.child_dict = new_child_dict

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            self.id = ''.join(new_id)


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

    def crop(self, *, remove_het_atoms: bool = True) -> None:
        """Crop to antigen binding domain.

        This method mutates the CD1 object.

        Args:
            remove_het_atoms: remove het atoms from structure as well

        """
        if hasattr(self, 'CD1'):
            alpha_chain = self.child_dict[self.CD1]

            new_alpha_chain = MHCchain(alpha_chain.id)

            for residue in alpha_chain:
                if residue.id[1] in IMGT_MH1_ABD or (not remove_het_atoms and residue.id[0] != ' '):
                    new_alpha_chain.add(residue.copy())

            new_alpha_chain.analyse(alpha_chain.chain_type)

            del self[alpha_chain.id]
            self.add(new_alpha_chain)

        if hasattr(self, 'B2M'):
            del self[self.B2M]
            self.B2M = None

    def standardise_chain_names(self) -> None:
        """Standardise CD1 chain name to A and B2M chain name to B."""
        new_id = []
        new_child_dict = {}

        if hasattr(self, 'CD1'):
            new_child_dict['A'] = self.child_dict[self.CD1]
            self.CD1 = 'A'
            new_id.append('A')

        if hasattr(self, 'B2M'):
            new_child_dict['B'] = self.child_dict[self.B2M]
            self.B2M = 'B'
            new_id.append('B')

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)

            for chain_id, chain in new_child_dict.items():
                chain.id = chain_id

        self.child_dict = new_child_dict

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            self.id = ''.join(new_id)


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

    def crop(self, *, remove_het_atoms: bool = True) -> None:
        """Crop to antigen binding domain.

        This method mutates the MR1 object.

        Args:
            remove_het_atoms: remove het atoms from structure as well

        """
        if hasattr(self, 'MR1'):
            alpha_chain = self.child_dict[self.CD1]

            new_alpha_chain = MHCchain(alpha_chain.id)

            for residue in alpha_chain:
                if residue.id[1] in IMGT_MH1_ABD or (not remove_het_atoms and residue.id[0] != ' '):
                    new_alpha_chain.add(residue.copy())

            new_alpha_chain.analyse(alpha_chain.chain_type)

            del self[alpha_chain.id]
            self.add(new_alpha_chain)

        if hasattr(self, 'B2M'):
            del self[self.B2M]
            self.B2M = None

    def standardise_chain_names(self) -> None:
        """Standardise MR1 chain name to A and B2M chain name to B."""
        new_id = []
        new_child_dict = {}

        if hasattr(self, 'MR1'):
            new_child_dict['A'] = self.child_dict[self.MR1]
            self.MR1 = 'A'
            new_id.append('A')

        if hasattr(self, 'B2M'):
            new_child_dict['B'] = self.child_dict[self.B2M]
            self.B2M = 'B'
            new_id.append('B')

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)

            for chain_id, chain in new_child_dict.items():
                chain.id = chain_id

        self.child_dict = new_child_dict

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            self.id = ''.join(new_id)


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

    def crop(self, *, remove_het_atoms: bool = True) -> None:
        """Crop to antigen binding domain.

        This method mutates the scMH1 object.

        Args:
            remove_het_atoms: remove het atoms from structure as well

        """
        alpha_chain = self.get_alpha()
        new_alpha_chain = MHCchain(alpha_chain.id)

        for residue in alpha_chain:
            if residue.id[1] in IMGT_MH1_ABD or (not remove_het_atoms and residue.id[0] != ' '):
                new_alpha_chain.add(residue.copy())

        new_alpha_chain.analyse(alpha_chain.chain_type)

        del self[alpha_chain.id]
        self.add(new_alpha_chain)

    def standardise_chain_names(self) -> None:
        """Standardise MHC chain name to A."""
        new_child_dict = {}
        for MH1_domain in set(['MH1', 'GA1', 'GA2']):
            if hasattr(self, MH1_domain):
                new_child_dict['A'] = self.child_dict[getattr(self, MH1_domain)]
                setattr(self, MH1_domain, 'A')
                break

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)

            for chain_id, chain in new_child_dict.items():
                chain.id = chain_id

        self.child_dict = new_child_dict

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            self.id = 'A'


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

    def crop(self, *, remove_het_atoms: bool = True) -> None:
        """Crop to antigen binding domain.

        This method mutates the CD1 object.

        Args:
            remove_het_atoms: remove het atoms from structure as well

        """
        if hasattr(self, 'CD1'):
            alpha_chain = self.child_dict[self.CD1]

            new_alpha_chain = MHCchain(alpha_chain.id)

            for residue in alpha_chain:
                if residue.id[1] in IMGT_MH1_ABD or (not remove_het_atoms and residue.id[0] != ' '):
                    new_alpha_chain.add(residue.copy())

            new_alpha_chain.analyse(alpha_chain.chain_type)

            del self[alpha_chain.id]
            self.add(new_alpha_chain)

    def standardise_chain_names(self) -> None:
        """Standardise CD1 chain name to A."""
        new_child_dict = {}
        if hasattr(self, 'CD1'):
            new_child_dict['A'] = self.child_dict[self.CD1]
            self.CD1 = 'A'

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)

            for chain_id, chain in new_child_dict.items():
                chain.id = chain_id

        self.child_dict = new_child_dict

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            self.id = 'A'


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

    def crop(self, *, remove_het_atoms: bool = True) -> None:
        """Crop to antigen binding domain.

        This method mutates the scMH2 object.

        Args:
            remove_het_atoms: remove het atoms from structure as well

        """
        new_child_dict = {}

        for chain in self:
            new_chain = MHCchain(chain.id)

            for residue in chain:
                if residue.id[1] in IMGT_MH2_ABD or (not remove_het_atoms and residue.id[0] != ' '):
                    new_chain.add(residue.copy())

            new_chain.analyse(chain.chain_type)
            new_child_dict[new_chain.id] = new_chain

        for chain_id in new_child_dict:
            del self[chain_id]

        for new_chain in new_child_dict.values():
            self.add(new_chain)

    def standardise_chain_names(self) -> None:
        """Standardise MHC chain 1 name to A or MHC chain 2 name to B."""
        new_id = []
        new_child_dict = {}

        if hasattr(self, 'GA'):
            new_child_dict['A'] = self.child_dict[self.GA]
            self.GA = 'A'
            new_id.append('A')

        if hasattr(self, 'GB'):
            new_child_dict['B'] = self.child_dict[self.GB]
            self.GB = 'B'
            new_id.append('B')

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)

            for chain_id, chain in new_child_dict.items():
                chain.id = chain_id

        self.child_dict = new_child_dict

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            self.id = ''.join(new_id)
