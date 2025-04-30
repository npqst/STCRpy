"""
constants.py
@author: leem
@date:   9 May 2017

Constant values that are useful. Based off of rotlib.constants.
"""

import re


def tuplefy(x):
    """
    Interpretation for converting numbering (in string) into a tuple.
    Args:
        x: A string for the identifier of a numbered position. e.g "H100A".
    Returns:
        A tuple of the chain tupe followed by a tuple of residue id and insertion code. eg. ( H, (100, "A") )

    """
    chain, resi, ins = re.split(r"(\d+)", x)
    assert chain in ["B", "A", "G", "D"], "Not a recognised chain type."
    ins = ins if ins else " "
    return (chain, (int(resi), ins))


#######
# String constants
#######

# Basic string names for residues
RESIDUES = {}
RESIDUES["ALA"] = "A"
RESIDUES["CYS"] = "C"
RESIDUES["ASP"] = "D"
RESIDUES["GLU"] = "E"
RESIDUES["PHE"] = "F"
RESIDUES["GLY"] = "G"
RESIDUES["HIS"] = "H"
RESIDUES["ILE"] = "I"
RESIDUES["LYS"] = "K"
RESIDUES["LEU"] = "L"
RESIDUES["MET"] = "M"
RESIDUES["ASN"] = "N"
RESIDUES["PRO"] = "P"
RESIDUES["GLN"] = "Q"
RESIDUES["ARG"] = "R"
RESIDUES["SER"] = "S"
RESIDUES["THR"] = "T"
RESIDUES["VAL"] = "V"
RESIDUES["TRP"] = "W"
RESIDUES["TYR"] = "Y"
RESIDUES_SINGLE = dict([(v, k) for k, v in list(RESIDUES.items())])
RESIDUES_SINGLE_STRING = "ACDEFGHIKLMNPQRSTVWY"

# Atoms in the backbone/CB for doing checks on residues.
BACKBONE_ATOMS = ["N", "CA", "C", "O"]
BACKBONE_CB = ["N", "CA", "C", "O", "CB"]

# TCR chain types
TCR_CHAINS = ["B", "A", "G", "D"]
TCR_REGIONS = [
    "fwb1",
    "fwb2",
    "fwb3",
    "fwb4",
    "fwa1",
    "fwa2",
    "fwa3",
    "fwa4",
    "fwg1",
    "fwg2",
    "fwg3",
    "fwg4",
    "fwd1",
    "fwd2",
    "fwd3",
    "fwd4",
    "cdrb1",
    "cdrb2",
    "cdrb3",
    "cdrb4",
    "cdra1",
    "cdra2",
    "cdra3",
    "cdra4",
    "cdrg1",
    "cdrg2",
    "cdrg3",
    "cdrg4",
    "cdrd1",
    "cdrd2",
    "cdrd3",
    "cdrd4",
]

# MHC chain types and regions
MHC_CHAINS = ["GA1", "GA2", "GA", "GB"]  # G-alpha1,2, G-alpha,beta
MHC_REGIONS = ["Astrand", "Bstrand", "Cstrand", "Dstrand", "Helix", "Turn"]

# Common names for species
COMMON_NAMES = {
    "bos taurus": "cattle",
    "camelus dromedarius": "arabian camel",
    "canis lupus familiaris": "domestic dog",
    "cerocebus atys": "sooty mangabey",
    "danio rerio": "zebrafish",
    "homo sapiens": "human",
    "macaca fascicularis": "crab-eating macaque",
    "macaca mulatta": "rhesus macaque",
    "macaca nemestrina": "Southern pig-tailed macaque",
    "mus musculus": "house mouse",
    "mus cookii": "cook's mouse",
    "mus minutoides": "African pygmy mouse",
    "mus pahari": "Gairdner's shrewmouse",
    "mus saxicola": "brown spiny mouse",
    "mus spretus": "Algerian mouse",
    "oncorhynchus mykiss": "rainbow trout",
    "ornithorhynchus anatinus": "platypus",
    "oryctolagus cuniculus": "rabbit",
    "ovis aries": "sheep",
    "papio anubis": "olive baboon",
    "staphylococcus aureus": "S. aureus",
    "rattus norvegicus": "norway rat",
    "rattus rattus": "black rat",
    "sus scrofa": "pig",
    "vicugna pacos": "alpaca",
}

# Atoms in the backbone/CB for doing checks on residues.
BACKBONE_ATOMS = ["N", "CA", "C", "O"]
BACKBONE_CB = ["N", "CA", "C", "O", "CB"]

SIDECHAIN_ATOMS = dict(
    ALA=["CB"],
    ARG=["CB", "CG", "CD", "CZ", "NE", "NH1", "NH2"],
    ASN=["CB", "CG", "ND2", "OD1"],
    ASP=["CB", "CG", "OD1", "OD2"],
    CYS=["CB", "SG"],
    GLN=["CB", "CG", "CD", "NE2", "OE1"],
    GLU=["CB", "CG", "CD", "OE1", "OE2"],
    HIS=["CB", "CG", "CD2", "CE1", "ND1", "NE2"],
    ILE=["CB", "CG1", "CD1", "CG2"],
    LEU=["CB", "CG", "CD1", "CD2"],
    LYS=["CB", "CG", "CD", "CE", "NZ"],
    MET=["CB", "CG", "SD", "CE"],
    PHE=["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    PRO=["CB", "CG", "CD"],
    SER=["CB", "OG"],
    THR=["CB", "CG2", "OG1"],
    TRP=["CB", "CG", "CD1", "CD2", "CE2", "CE3", "CH2", "CZ2", "CZ3", "NE1"],
    TYR=["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    VAL=["CB", "CG1", "CG2"],
)

# BLOSUM62 matrix
BLOSUM62 = {
    ("A", "A"): 4,
    ("C", "A"): 0,
    ("C", "C"): 9,
    ("C", "D"): -3,
    ("C", "N"): -3,
    ("C", "R"): -3,
    ("D", "A"): -2,
    ("D", "D"): 6,
    ("D", "N"): 1,
    ("D", "R"): -2,
    ("E", "A"): -1,
    ("E", "C"): -4,
    ("E", "D"): 2,
    ("E", "E"): 5,
    ("E", "N"): 0,
    ("E", "Q"): 2,
    ("E", "R"): 0,
    ("F", "A"): -2,
    ("F", "C"): -2,
    ("F", "D"): -3,
    ("F", "E"): -3,
    ("F", "F"): 6,
    ("F", "G"): -3,
    ("F", "H"): -1,
    ("F", "I"): 0,
    ("F", "K"): -3,
    ("F", "L"): 0,
    ("F", "M"): 0,
    ("F", "N"): -3,
    ("F", "Q"): -3,
    ("F", "R"): -3,
    ("G", "A"): 0,
    ("G", "C"): -3,
    ("G", "D"): -1,
    ("G", "E"): -2,
    ("G", "G"): 6,
    ("G", "N"): 0,
    ("G", "Q"): -2,
    ("G", "R"): -2,
    ("H", "A"): -2,
    ("H", "C"): -3,
    ("H", "D"): -1,
    ("H", "E"): 0,
    ("H", "G"): -2,
    ("H", "H"): 8,
    ("H", "N"): 1,
    ("H", "Q"): 0,
    ("H", "R"): 0,
    ("I", "A"): -1,
    ("I", "C"): -1,
    ("I", "D"): -3,
    ("I", "E"): -3,
    ("I", "G"): -4,
    ("I", "H"): -3,
    ("I", "I"): 4,
    ("I", "N"): -3,
    ("I", "Q"): -3,
    ("I", "R"): -3,
    ("K", "A"): -1,
    ("K", "C"): -3,
    ("K", "D"): -1,
    ("K", "E"): 1,
    ("K", "G"): -2,
    ("K", "H"): -1,
    ("K", "I"): -3,
    ("K", "K"): 5,
    ("K", "L"): -2,
    ("K", "N"): 0,
    ("K", "Q"): 1,
    ("K", "R"): 2,
    ("L", "A"): -1,
    ("L", "C"): -1,
    ("L", "D"): -4,
    ("L", "E"): -3,
    ("L", "G"): -4,
    ("L", "H"): -3,
    ("L", "I"): 2,
    ("L", "L"): 4,
    ("L", "N"): -3,
    ("L", "Q"): -2,
    ("L", "R"): -2,
    ("M", "A"): -1,
    ("M", "C"): -1,
    ("M", "D"): -3,
    ("M", "E"): -2,
    ("M", "G"): -3,
    ("M", "H"): -2,
    ("M", "I"): 1,
    ("M", "K"): -1,
    ("M", "L"): 2,
    ("M", "M"): 5,
    ("M", "N"): -2,
    ("M", "Q"): 0,
    ("M", "R"): -1,
    ("N", "A"): -2,
    ("N", "N"): 6,
    ("N", "R"): 0,
    ("P", "A"): -1,
    ("P", "C"): -3,
    ("P", "D"): -1,
    ("P", "E"): -1,
    ("P", "F"): -4,
    ("P", "G"): -2,
    ("P", "H"): -2,
    ("P", "I"): -3,
    ("P", "K"): -1,
    ("P", "L"): -3,
    ("P", "M"): -2,
    ("P", "N"): -2,
    ("P", "P"): 7,
    ("P", "Q"): -1,
    ("P", "R"): -2,
    ("Q", "A"): -1,
    ("Q", "C"): -3,
    ("Q", "D"): 0,
    ("Q", "N"): 0,
    ("Q", "Q"): 5,
    ("Q", "R"): 1,
    ("R", "A"): -1,
    ("R", "R"): 5,
    ("S", "A"): 1,
    ("S", "C"): -1,
    ("S", "D"): 0,
    ("S", "E"): 0,
    ("S", "F"): -2,
    ("S", "G"): 0,
    ("S", "H"): -1,
    ("S", "I"): -2,
    ("S", "K"): 0,
    ("S", "L"): -2,
    ("S", "M"): -1,
    ("S", "N"): 1,
    ("S", "P"): -1,
    ("S", "Q"): 0,
    ("S", "R"): -1,
    ("S", "S"): 4,
    ("T", "A"): 0,
    ("T", "C"): -1,
    ("T", "D"): -1,
    ("T", "E"): -1,
    ("T", "F"): -2,
    ("T", "G"): -2,
    ("T", "H"): -2,
    ("T", "I"): -1,
    ("T", "K"): -1,
    ("T", "L"): -1,
    ("T", "M"): -1,
    ("T", "N"): 0,
    ("T", "P"): -1,
    ("T", "Q"): -1,
    ("T", "R"): -1,
    ("T", "S"): 1,
    ("T", "T"): 5,
    ("V", "A"): 0,
    ("V", "C"): -1,
    ("V", "D"): -3,
    ("V", "E"): -2,
    ("V", "F"): -1,
    ("V", "G"): -3,
    ("V", "H"): -3,
    ("V", "I"): 3,
    ("V", "K"): -2,
    ("V", "L"): 1,
    ("V", "M"): 1,
    ("V", "N"): -3,
    ("V", "P"): -2,
    ("V", "Q"): -2,
    ("V", "R"): -3,
    ("V", "S"): -2,
    ("V", "T"): 0,
    ("V", "V"): 4,
    ("V", "W"): -3,
    ("V", "Y"): -1,
    ("W", "A"): -3,
    ("W", "C"): -2,
    ("W", "D"): -4,
    ("W", "E"): -3,
    ("W", "F"): 1,
    ("W", "G"): -2,
    ("W", "H"): -2,
    ("W", "I"): -3,
    ("W", "K"): -3,
    ("W", "L"): -2,
    ("W", "M"): -1,
    ("W", "N"): -4,
    ("W", "P"): -4,
    ("W", "Q"): -2,
    ("W", "R"): -3,
    ("W", "S"): -3,
    ("W", "T"): -2,
    ("W", "W"): 11,
    ("Y", "A"): -2,
    ("Y", "C"): -2,
    ("Y", "D"): -3,
    ("Y", "E"): -2,
    ("Y", "F"): 3,
    ("Y", "G"): -3,
    ("Y", "H"): 2,
    ("Y", "I"): -1,
    ("Y", "K"): -2,
    ("Y", "L"): -1,
    ("Y", "M"): -1,
    ("Y", "N"): -2,
    ("Y", "P"): -3,
    ("Y", "Q"): -1,
    ("Y", "R"): -2,
    ("Y", "S"): -2,
    ("Y", "T"): -2,
    ("Y", "W"): 2,
    ("Y", "Y"): 7,
}
