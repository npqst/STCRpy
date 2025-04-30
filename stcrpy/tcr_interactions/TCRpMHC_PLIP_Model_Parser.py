import os
import warnings
import copy

try:
    from plip.structure.preparation import PDBComplex
except ModuleNotFoundError:
    warnings.warn(
        """\n\nPLIP package not found. \nProfiling interactions will not be possible \nTo enable interaction profiling, install PLIP with:
        \npip install plip --no-deps\n\n"""
    )

from rdkit import Chem
from Bio import PDB
from Bio.PDB.PDBParser import PDBParser

from ..tcr_processing.TCRParser import TCRParser
from ..tcr_processing.TCR import TCR


class TCRpMHC_PLIP_Model_Parser:
    def __init__(self, tmp_dir=None):
        self.parser = PDBParser()
        self.tcr_parser = TCRParser()
        self.tmp_dir = tmp_dir if tmp_dir is not None else "./"

    def parse_tcr_pmhc_complex(
        self,
        tcr_pmhc_complex: TCR,
        delete_tmp_files=True,
        renumber=True,
    ) -> "PDBComplex":

        # tcr_pmhc_complex = copy.deepcopy(
        #     tcr_pmhc_complex
        # )  # copy the complex to prevent renumbering from persisting in TCR object

        ligand = PDB.Model.Model(id=0)

        peptide_chain = tcr_pmhc_complex.antigen
        if len(peptide_chain) != 1:
            # try to identify correct antigen
            # if a chain is longer than 25 residues reject it, this may happen if anarci has failed to label an MHC chain
            peptide_chain = [
                c
                for c in peptide_chain
                if (
                    len(c) < 25
                    and isinstance(c, PDB.Chain.Chain)
                    or isinstance(c, PDB.Residue.Residue)
                )
            ]

        assert (
            len(peptide_chain) == 1
        ), f"More or less than one peptide chain found: {peptide_chain}"

        if isinstance(
            peptide_chain[0], PDB.Residue.Residue
        ):  # wrap single residue antigen in chain
            residue_as_chain = PDB.Chain.Chain(id="Z")
            residue_as_chain.add(peptide_chain[0].copy())
            peptide_chain = [residue_as_chain]

        ligand.add(peptide_chain[0].copy())

        tcr_and_mhc_chains = [
            c.copy()
            for c in list(tcr_pmhc_complex.get_chains())
            + list(tcr_pmhc_complex.get_MHC()[0].get_chains())
        ]
        if renumber:
            # renumber each chain from one to N to avoid automated renumbering issues related to plip and openbabel
            renumbering = {}
            for chain in tcr_and_mhc_chains:
                renumbering[chain.id] = {}
                for new_idx, res in enumerate(chain.get_residues()):
                    new_id = (" ", new_idx + 1, " ")
                    renumbering[chain.id][new_id] = res.id
                    res.id = new_id
            domain_assignment = tcr_pmhc_complex.get_domain_assignment()

        TCR_MHC_FILE = os.path.join(self.tmp_dir, "tcr_mhc.pdb")
        PEPTIDE_PDB_FILE = os.path.join(self.tmp_dir, "peptide.pdb")
        PEPTIDE_SDF_FILE = os.path.join(self.tmp_dir, "peptide.sdf")

        io = PDB.PDBIO()
        io.set_structure(ligand)
        io.save(PEPTIDE_PDB_FILE)

        tcr_mhc_struct = PDB.Model.Model(id=0)
        # add TCR chains to protein structure
        for chain in tcr_pmhc_complex.get_chains():
            tcr_mhc_struct.add(chain)
        # add MHC chain to protein structure
        for chain in tcr_pmhc_complex.get_MHC()[0].get_chains():
            tcr_mhc_struct.add(chain)

        io = PDB.PDBIO()
        io.set_structure(tcr_mhc_struct)
        io.save(TCR_MHC_FILE)

        rdkit_mol = Chem.MolFromPDBFile(PEPTIDE_PDB_FILE)
        Chem.MolToMolFile(rdkit_mol, PEPTIDE_SDF_FILE)
        with open(TCR_MHC_FILE, "r") as f:
            protein = f.read()
        protein = [line for line in protein.split("\n") if line.startswith("ATOM")]
        ligand = Chem.MolFromMolFile(PEPTIDE_SDF_FILE)
        ligand_pdb_block = Chem.MolToPDBBlock(ligand)
        complex_pdb_block = "\n".join(protein) + "\n" + ligand_pdb_block
        complex = PDBComplex()
        complex.load_pdb(complex_pdb_block, as_string=True)

        if delete_tmp_files:
            os.remove(TCR_MHC_FILE)
            os.remove(PEPTIDE_PDB_FILE)
            os.remove(PEPTIDE_SDF_FILE)

            if renumber:
                return complex, renumbering, domain_assignment
            else:
                return complex

        else:
            if renumber:
                return (
                    complex,
                    TCR_MHC_FILE,
                    PEPTIDE_PDB_FILE,
                    PEPTIDE_SDF_FILE,
                    renumbering,
                    domain_assignment,
                )
            else:
                return complex, TCR_MHC_FILE, PEPTIDE_PDB_FILE, PEPTIDE_SDF_FILE

    def map_amino_acids_to_ligands(self, ligand_pdb, ligand_sdf):
        ligand_structure = self.parser.get_structure("tmp", ligand_pdb)
        sdf_supplier = Chem.SDMolSupplier(ligand_sdf)
        mol = [x for x in sdf_supplier][0]
        sdf_coords = mol.GetConformer().GetPositions()
        coord_to_aa = {}
        for coord in sdf_coords:
            pdb_atom = [
                a
                for a in ligand_structure.get_atoms()
                if sum((a.coord - coord) ** 2) < 0.0001
            ]
            assert len(pdb_atom) == 1
            coord_to_aa[tuple(coord)] = (
                pdb_atom[0].parent.resname,
                pdb_atom[0].parent.id[1],
            )
        return coord_to_aa
