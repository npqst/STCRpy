
import os
import numpy as np
import Bio

from ..tcr_processing.TCRParser import TCRParser
from ..tcr_processing.TCRIO import TCRIO
from ..tcr_processing import TCR, MHCchain


class TCRCoM():
    def __init__(self):
        self.set_reffile()
        
        tcr_parser = TCRParser()
        self.ref_model = list(tcr_parser.get_tcr_structure("reference", self.reffile).get_TCRs())[0]
        self.set_reference_residues()
    
    def set_reffile(self):
        raise NotImplementedError('TCRCom cannot be insantiated directly, instantiate its subclass')

    def set_reference_residues(self):
        self.set_tcr_reference()
        self.set_mhc_reference()

    def set_tcr_reference(self):
        self.reference_VA_residues = [r for r in self.ref_model.get_VA().get_residues() if r.get_id()[1] >= 1 and r.get_id()[1] <= 121]
        self.reference_VB_residues = [r for r in self.ref_model.get_VB().get_residues() if r.get_id()[1] >= 1 and r.get_id()[1] <= 126]
        
    def set_mhc_reference(self):
        raise NotImplementedError('TCRCom cannot be insantiated directly, instantiate its subclass')

    def get_filtered_TCR_residues(self):
        tcr_A_residues = [self.tcr.get_VA()[r.get_id()] for r in self.reference_VA_residues if r.get_id() in self.tcr.get_VA()]
        tcr_B_residues = [self.tcr.get_VB()[r.get_id()] for r in self.reference_VB_residues if r.get_id() in self.tcr.get_VB()]
        return tcr_A_residues + tcr_B_residues

    def center_of_mass(self, entity, geometric=False):
        # Structure, Model, Chain, Residue
        if isinstance(entity, Bio.PDB.Entity.Entity):
            atom_list = entity.get_atoms()
        # List of Atoms
        elif hasattr(entity, "__iter__") and [x for x in entity if x.level == "A"]:
            atom_list = entity
        else:
            raise ValueError(
                "Center of Mass can only be calculated from the following objects:\n"
                "Structure, Model, Chain, Residue, list of Atoms."
            )
        
        masses, positions = zip(*[(atom.mass, atom.coord) for atom in atom_list])
        positions = np.asarray(positions)
        if "ukn" in set(masses) and not geometric:
            raise ValueError(
                "Some Atoms don't have an element assigned.\n"
                "Try adding them manually or calculate the geometrical center of mass instead"
            )
        
        if geometric:
            return positions.sum(axis=0) / len(atom_list)
        else:
            return np.matmul(np.asarray(masses), positions) / len(atom_list)

    def add_com(self, mhc_com, tcr_com):
        """
        Function to add pseudoatoms at MHC-CoM, TCR-CoM, and XYZ axis to the output PDB file
        """
        new_structure = self.tcr.copy()
        
        # mhc_com
        mhc_com_chain = "X"
        new_structure.add(Bio.PDB.Chain.Chain(mhc_com_chain))
        res_id = (" ", 1, " ")
        new_residue = Bio.PDB.Residue.Residue(res_id, "MCM", " ")
        new_atom = Bio.PDB.Atom.Atom("C", mhc_com, 0, 0.0, " ", "C", 1, "C")
        new_residue.add(new_atom)
        new_structure.child_dict[mhc_com_chain].add(new_residue)
        
        # tcr com
        tcr_com_chain = "Y"
        new_structure.add(Bio.PDB.Chain.Chain(tcr_com_chain))
        res_id = (" ", 1, " ")
        new_residue = Bio.PDB.Residue.Residue(res_id, "TCM", " ")
        new_atom = Bio.PDB.Atom.Atom("C", tcr_com, 0, 0.0, " ", "C", 1, "C")
        new_residue.add(new_atom)
        new_structure.child_dict[tcr_com_chain].add(new_residue)
        
        # X,Y,Z atoms
        pos = [[50, 0, 0], [0, 50, 0], [0, 0, 50]]
        resn = ["X", "Y", "Z"]
        xyz_chain = "Z"
        new_structure.add(Bio.PDB.Chain.Chain(xyz_chain))
        for i in [0, 1, 2]:
            res_id = (" ", i + 1, " ")
            new_residue = Bio.PDB.Residue.Residue(res_id, resn[i], " ")
            new_atom = Bio.PDB.Atom.Atom("O", pos[i], 0, 0.0, " ", "O", 1, "O")
            new_residue.add(new_atom)
            new_structure.child_dict[xyz_chain].add(new_residue)
        
        return new_structure

    def calculate_geometry(
            self,
            tcr: TCR,
            save_aligned_as: str = None,
    ):
        self.tcr = tcr

        ref_mhc_atoms = [res["CA"] for res in self.reference_MHC_residues]

        mhc_atoms = [res["CA"] for res in self.get_filtered_MHC_residues()]

        superimposer = Bio.PDB.Superimposer()
        superimposer.set_atoms(ref_mhc_atoms, mhc_atoms)
        superimposer.apply(self.tcr.parent.get_atoms())

        mhc_com = self.center_of_mass(mhc_atoms, geometric=True)

        tcr_atoms = [res["CA"] for res in self.get_filtered_TCR_residues()]
        tcr_com = self.center_of_mass(tcr_atoms, geometric=True)

        if save_aligned_as:
            aligned_tcr = self.add_com(mhc_com, tcr_com)
            io = TCRIO()
            io.save(aligned_tcr, save_as=save_aligned_as)
        
        com_distance = tcr_com - mhc_com
        self.r = np.sqrt(np.sum(np.square(com_distance)))
        self.theta = np.degrees(np.arctan2(com_distance[1], com_distance[0]))
        self.phi = np.degrees(np.arccos(com_distance[2] / self.r))

        return self.r, self.theta, self.phi


class MHCI_TCRCoM(TCRCoM):
    def __init__(self):
        super().__init__()
    
    def set_reffile(self):
        """Sets reference file for MHC class I structures
        """
        self.reffile = os.path.join(os.path.dirname(__file__), 'reference_data/dock_reference_1_imgt_numbered.pdb')

    def set_mhc_reference(self):
        mhc = self.ref_model.get_MHC()[0].get_MH1()
        self.reference_MHC_residues = [r for r in mhc.get_residues() if r.get_id()[1] >= 1 and r.get_id()[1] <= 179]

    def get_filtered_MHC_residues(self):
        mhc = self.tcr.get_MHC()[0]
        if not isinstance(mhc, MHCchain):           # handle single MHC chain case
            mhc = mhc.get_MH1()
        filtered_MHC_residues = [mhc[ref_res.get_id()] for ref_res in self.reference_MHC_residues if ref_res.get_id() in mhc]
        return filtered_MHC_residues
    

class MHCII_TCRCoM(TCRCoM):
    def __init__(self):
        super().__init__()
    
    def set_reffile(self):
        """Sets reference file for MHC class II structures
        """
        self.reffile = os.path.join(os.path.dirname(__file__), 'reference_data/dock_reference_2_imgt_numbered.pdb')

    def set_mhc_reference(self):
        mhc = self.ref_model.get_MHC()[0]
        self.reference_MHC_residues = [r for r in mhc.get_GA().get_residues() if r.get_id()[1] >= 1 and r.get_id()[1] <= 88]
        self.reference_MHC_residues.extend([r for r in mhc.get_GB().get_residues() if r.get_id()[1] >= 1 and r.get_id()[1] <= 87])
    
    def get_filtered_MHC_residues(self):
        mhc = self.tcr.get_MHC()[0]
        filtered_MHC_residues = [
            mhc.get_GA()[ref_res.get_id()] 
            for ref_res in self.reference_MHC_residues 
            if ref_res.parent.chain_type == 'GA' and ref_res.get_id() in mhc.get_GA()
            ]
        filtered_MHC_residues.extend([
            mhc.get_GB()[ref_res.get_id()] 
            for ref_res in self.reference_MHC_residues
            if ref_res.parent.chain_type == 'GB' and ref_res.get_id() in mhc.get_GB()
            ])
        return filtered_MHC_residues
