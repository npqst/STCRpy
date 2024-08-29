import warnings

try:
    from plip.structure.preparation import PDBComplex
except ModuleNotFoundError:
    warnings.warn(
        """\n\nPLIP package not found. \nProfiling interactions will not be possible \nTo enable interaction profiling, install PLIP with:
        \npip install plip --no-deps\n\n"""
    )
from rdkit import Chem


def return_interactions(
    protein_file=None, ligand_file=None, complex_file=None, pymol_visualization=False
):
    with open(protein_file, "r") as f:
        protein = f.read()
    protein = [line for line in protein.split("\n") if line.startswith("ATOM")]
    ligand = Chem.MolFromMolFile(ligand_file)
    ligand_pdb_block = Chem.MolToPDBBlock(ligand)
    complex_pdb_block = "\n".join(protein) + "\n" + ligand_pdb_block
    # return complex_pdb_block, ligand_pdb_block, protein
    my_mol = PDBComplex()
    my_mol.load_pdb(complex_pdb_block, as_string=True)
    my_mol.analyze()
    return my_mol


class Interaction:

    def __init__(
        self,
        type,
        protein_atom,
        protein_chain,
        protein_residue,
        protein_number,
        ligand_atom,
        distance,
        angle,
        plip_id,
    ) -> None:
        self.type = type
        self.protein_atom = protein_atom
        self.protein_chain = protein_chain
        self.protein_residue = protein_residue
        self.protein_number = protein_number
        self.ligand_atom = ligand_atom
        self.distance = distance
        self.angle = angle
        self.plip_id = plip_id

    def to_tuple(self):
        return (
            self.type,
            self.protein_atom,
            self.protein_chain,
            self.protein_residue,
            self.protein_number,
            self.ligand_atom,
            self.distance,
            self.angle,
            self.plip_id,
        )


def parse_interaction(interaction) -> Interaction:
    if "saltbridge" in str(type(interaction)):
        return Interaction("saltbridge", *process_saltbridge(interaction))
    elif "hydroph" in str(type(interaction)):
        return Interaction("hydrophobic", *process_hydrophobic(interaction))
    elif "hbond" in str(type(interaction)):
        return Interaction("hbond", *process_hbond(interaction))
    elif "pistack" in str(type(interaction)):
        return Interaction("pistack", *process_pi_stack(interaction))
    else:
        raise NotImplementedError(f"Parsing not implemented for {type(interaction)}")


def process_pi_stack(interaction):
    protein_ring_atoms = [
        (j.coords, j.atomicnum) for j in interaction.proteinring.atoms
    ]
    protein_chain = interaction.reschain
    protein_residue = interaction.restype
    protein_number = interaction.resnr
    ligand_ring_atoms = [(j.coords, j.atomicnum) for j in interaction.ligandring.atoms]
    distance = interaction.distance
    angle = interaction.angle
    plip_id = None
    return (
        protein_ring_atoms,
        protein_chain,
        protein_residue,
        protein_number,
        ligand_ring_atoms,
        distance,
        angle,
        plip_id,
    )


def process_hydrophobic(interaction):
    protein_atom = [(interaction.bsatom.coords, interaction.bsatom.atomicnum)]
    protein_chain = interaction.reschain
    protein_residue = interaction.restype
    protein_number = interaction.resnr
    ligand_atom = [(interaction.ligatom.coords, interaction.ligatom.atomicnum)]
    distance = interaction.distance
    plip_id = None
    return (
        protein_atom,
        protein_chain,
        protein_residue,
        protein_number,
        ligand_atom,
        distance,
        None,
        plip_id,
    )


def process_hbond(interaction):
    if interaction.protisdon:
        protein_atom = [(interaction.d.coords, interaction.d.atomicnum)]
        ligand_atom = [(interaction.a.coords, interaction.a.atomicnum)]
    else:
        protein_atom = [(interaction.a.coords, interaction.a.atomicnum)]
        ligand_atom = [(interaction.d.coords, interaction.d.atomicnum)]

    protein_chain = interaction.reschain
    protein_residue = interaction.restype
    protein_number = interaction.resnr
    distance = interaction.distance_ad
    angle = interaction.angle
    plip_id = None
    return (
        protein_atom,
        protein_chain,
        protein_residue,
        protein_number,
        ligand_atom,
        distance,
        angle,
        plip_id,
    )


def process_saltbridge(interaction):
    if interaction.protispos:
        protein_atom = [(a.coords, a.atomicnum) for a in interaction.positive.atoms]
        ligand_atom = [(a.coords, a.atomicnum) for a in interaction.negative.atoms]
    else:
        protein_atom = [(a.coords, a.atomicnum) for a in interaction.negative.atoms]
        ligand_atom = [(a.coords, a.atomicnum) for a in interaction.positive.atoms]
    protein_chain = interaction.reschain
    protein_residue = interaction.restype
    protein_number = interaction.resnr
    distance = interaction.distance
    plip_id = None
    return (
        protein_atom,
        protein_chain,
        protein_residue,
        protein_number,
        ligand_atom,
        distance,
        None,
        plip_id,
    )
