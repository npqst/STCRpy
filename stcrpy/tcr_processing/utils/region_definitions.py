"""
A module to deal with region annotations for IMGT scheme.
"""

from .constants import TCR_CHAINS, TCR_REGIONS

IMGT_CDR_BOUNDARIES = {
    "1": {"imgt": (27, 38)},
    "2": {"imgt": (56, 65)},
    "3": {"imgt": (105, 117)},
}

# regions for TCR
_regions = {"imgt": {}}
_regions["imgt"]["A"] = _regions["imgt"]["B"] = (
    "11111111111111111111111111222222222222333333333333333334444444444555555555555555555555555555555555555555666666666666677777777777"
)

# Set the IMGT definitions for TCR chain types
_regions["imgt"]["G"] = _regions["imgt"]["B"]
_regions["imgt"]["D"] = _regions["imgt"]["A"]


# For internal use only. These are not direct conversions and are handled heuristically.
# Currently only using the IMGT numbering scheme although other numbering schemes may be introduced.
_index_to_imgt_state = {
    ("imgt", "B"): {
        1: 0,
        2: 1,
        3: 2,
        4: 3,
        5: 4,
        6: 5,
        7: 6,
        8: 7,
        9: 8,
        10: 9,
        11: 10,
        12: 11,
        13: 12,
        14: 13,
        15: 14,
        16: 15,
        17: 16,
        18: 17,
        19: 18,
        20: 19,
        21: 20,
        22: 21,
        23: 22,
        24: 23,
        25: 24,
        26: 25,
        27: 26,
        28: 27,
        29: 28,
        30: 29,
        31: 30,
        32: 31,
        33: 32,
        34: 33,
        35: 34,
        36: 35,
        37: 36,
        38: 37,
        39: 38,
        40: 39,
        41: 40,
        42: 41,
        43: 42,
        44: 43,
        45: 44,
        46: 45,
        47: 46,
        48: 47,
        49: 48,
        50: 49,
        51: 50,
        52: 51,
        53: 52,
        54: 53,
        55: 54,
        56: 55,
        57: 56,
        58: 57,
        59: 58,
        60: 59,
        61: 60,
        62: 61,
        63: 62,
        64: 63,
        65: 64,
        66: 65,
        67: 66,
        68: 67,
        69: 68,
        70: 69,
        71: 70,
        72: 71,
        73: 72,
        74: 73,
        75: 74,
        76: 75,
        77: 76,
        78: 77,
        79: 78,
        80: 79,
        81: 80,
        82: 81,
        83: 82,
        84: 83,
        85: 84,
        86: 85,
        87: 86,
        88: 87,
        89: 88,
        90: 89,
        91: 90,
        92: 91,
        93: 92,
        94: 93,
        95: 94,
        96: 95,
        97: 96,
        98: 97,
        99: 98,
        100: 99,
        101: 100,
        102: 101,
        103: 102,
        104: 103,
        105: 104,
        106: 105,
        107: 106,
        108: 107,
        109: 108,
        110: 109,
        111: 110,
        112: 111,
        113: 112,
        114: 113,
        115: 114,
        116: 115,
        117: 116,
        118: 117,
        119: 118,
        120: 119,
        121: 120,
        122: 121,
        123: 122,
        124: 123,
        125: 124,
        126: 125,
        127: 126,
        128: 127,
    },
    ("imgt", "A"): {
        1: 0,
        2: 1,
        3: 2,
        4: 3,
        5: 4,
        6: 5,
        7: 6,
        8: 7,
        9: 8,
        10: 9,
        11: 10,
        12: 11,
        13: 12,
        14: 13,
        15: 14,
        16: 15,
        17: 16,
        18: 17,
        19: 18,
        20: 19,
        21: 20,
        22: 21,
        23: 22,
        24: 23,
        25: 24,
        26: 25,
        27: 26,
        28: 27,
        29: 28,
        30: 29,
        31: 30,
        32: 31,
        33: 32,
        34: 33,
        35: 34,
        36: 35,
        37: 36,
        38: 37,
        39: 38,
        40: 39,
        41: 40,
        42: 41,
        43: 42,
        44: 43,
        45: 44,
        46: 45,
        47: 46,
        48: 47,
        49: 48,
        50: 49,
        51: 50,
        52: 51,
        53: 52,
        54: 53,
        55: 54,
        56: 55,
        57: 56,
        58: 57,
        59: 58,
        60: 59,
        61: 60,
        62: 61,
        63: 62,
        64: 63,
        65: 64,
        66: 65,
        67: 66,
        68: 67,
        69: 68,
        70: 69,
        71: 70,
        72: 71,
        73: 72,
        74: 73,
        75: 74,
        76: 75,
        77: 76,
        78: 77,
        79: 78,
        80: 79,
        81: 80,
        82: 81,
        83: 82,
        84: 83,
        85: 84,
        86: 85,
        87: 86,
        88: 87,
        89: 88,
        90: 89,
        91: 90,
        92: 91,
        93: 92,
        94: 93,
        95: 94,
        96: 95,
        97: 96,
        98: 97,
        99: 98,
        100: 99,
        101: 100,
        102: 101,
        103: 102,
        104: 103,
        105: 104,
        106: 105,
        107: 106,
        108: 107,
        109: 108,
        110: 109,
        111: 110,
        112: 111,
        113: 112,
        114: 113,
        115: 114,
        116: 115,
        117: 116,
        118: 117,
        119: 118,
        120: 119,
        121: 120,
        122: 121,
        123: 122,
        124: 123,
        125: 124,
        126: 125,
        127: 126,
        128: 127,
    },
}

# IMGT states are the same across the board for B/D and A/G
_index_to_imgt_state[("imgt", "G")] = _index_to_imgt_state[("imgt", "B")]
_index_to_imgt_state[("imgt", "D")] = _index_to_imgt_state[("imgt", "A")]
_reg_one2three = {
    "1": "fw%s1",
    "2": "cdr%s1",
    "3": "fw%s2",
    "4": "cdr%s2",
    "5": "fw%s3",
    "6": "cdr%s3",
    "7": "fw%s4",
}

# regions for MHC
# This is based on the state vector used for numbering MHCs in ANARCI
# http://www.imgt.org/IMGTrepertoireMHC/Proteins/protein/G-DOMAIN/Gdomains.html
# Refer to IMGT domain align tool as well
#           1     7   10 14   18 21    28  31    38   42 45 49       50  54        61      68  72  74    80        90
# 0987654321|.....|A..|...|...|..|......|..|......|...|..|...|1234567|...|ABC......|AB......|...|A.|.....|.........|..
mhc_svec = "1111111111111116662222222222266333333336664444444444444445555555555555555555555555555555555555555555566666"
# ^ This is the state vector for mhc g-domains for regions

# Manual observation for residues 88-92: this seems like a loopy region? keep it as "turn" for now.
# state vector for CD1/MR1 GA1-like and GA2-like domains
#           1     7   10 14   18 21    28  31    38   42 45 49       50  54      61      68  72  74    80        90
# 0987654321|.....|A..|...|...|..|......|..|......|...|..|...|1234567|...|A......|AB......|...|A.|.....|.........|..
cd1_svec = "11111111111111166622222222222663333333366644444444444444455555555555555555555555555555555555555555555555"
# ^ This is the state vector for CD1 ga-like-domains for regions


# This is the state vector for mhc c-domains (B2M)
#           1        10   15  16     232627      383941  45      77     84             85  89     96 97   104105       117118
#   87654321|........|....|123|......|..||........||.|...|1234567|......|12345677654321|...|......|12|......||...........||.........|
mhc_cvec = "1111111111111119992222222222299999999993333333999999944444444999999999999995555555555559966666666999999999999977777777777"

_reg_tostring = {
    "1": "Astrand",
    "2": "Bstrand",
    "3": "Cstrand",
    "4": "Dstrand",
    "5": "Helix",
    "6": "Turn",
}  # http://www.imgt.org/IMGTScientificChart/Numbering/IMGTGsuperfamily.html
_reg_mhc_cdom = {
    "1": "Astrand",
    "2": "Bstrand",
    "3": "Cstrand",
    "4": "Dstrand",
    "5": "Estrand",
    "6": "Fstrand",
    "7": "Gstrand",
    "9": "Turn",
}  # C-Like domain

_regions["imgt"]["MH1"] = mhc_svec
_regions["imgt"]["CD1"] = cd1_svec
_regions["imgt"]["MR1"] = cd1_svec
_regions["imgt"]["GA"] = mhc_svec
_regions["imgt"]["GB"] = _regions["imgt"]["GA"]
_regions["imgt"]["GA1"] = mhc_svec
_regions["imgt"]["GA2"] = mhc_svec
_regions["imgt"]["GA1L"] = cd1_svec
_regions["imgt"]["GA2L"] = cd1_svec

# C-LIKE And B2-Microglobulin regions
_regions["imgt"]["B2M"] = mhc_cvec

_index_to_imgt_state[("imgt", "MH1")] = {
    1: 0,
    2: 1,
    3: 2,
    4: 3,
    5: 4,
    6: 5,
    7: 6,
    8: 8,
    9: 9,
    10: 10,
    11: 11,
    12: 12,
    13: 13,
    14: 14,
    15: 15,
    16: 16,
    17: 17,
    18: 18,
    19: 19,
    20: 20,
    21: 21,
    22: 22,
    23: 23,
    24: 24,
    25: 25,
    26: 26,
    27: 27,
    28: 28,
    29: 29,
    30: 30,
    31: 31,
    32: 32,
    33: 33,
    34: 34,
    35: 35,
    36: 36,
    37: 37,
    38: 38,
    39: 39,
    40: 40,
    41: 41,
    42: 42,
    43: 43,
    44: 44,
    45: 45,
    46: 46,
    47: 47,
    48: 48,
    49: 49,
    50: 57,
    51: 58,
    52: 59,
    53: 60,
    54: 61,
    55: 65,
    56: 66,
    57: 67,
    58: 68,
    59: 69,
    60: 70,
    61: 71,
    62: 74,
    63: 75,
    64: 76,
    65: 77,
    66: 78,
    67: 79,
    68: 80,
    69: 81,
    70: 82,
    71: 83,
    72: 84,
    73: 86,
    74: 87,
    75: 88,
    76: 89,
    77: 90,
    78: 91,
    79: 92,
    80: 93,
    81: 94,
    82: 95,
    83: 96,
    84: 97,
    85: 98,
    86: 99,
    87: 100,
    88: 101,
    89: 102,
    90: 103,
    91: 104,
    92: 105,
}
_index_to_imgt_state[("imgt", "CD1")] = {
    1: 0,
    2: 1,
    3: 2,
    4: 3,
    5: 4,
    6: 5,
    7: 6,
    8: 8,
    9: 9,
    10: 10,
    11: 11,
    12: 12,
    13: 13,
    14: 14,
    15: 15,
    16: 16,
    17: 17,
    18: 18,
    19: 19,
    20: 20,
    21: 21,
    22: 22,
    23: 23,
    24: 24,
    25: 25,
    26: 26,
    27: 27,
    28: 28,
    29: 29,
    30: 30,
    31: 31,
    32: 32,
    33: 33,
    34: 34,
    35: 35,
    36: 36,
    37: 37,
    38: 38,
    39: 39,
    40: 40,
    41: 41,
    42: 42,
    43: 43,
    44: 44,
    45: 45,
    46: 46,
    47: 47,
    48: 48,
    49: 49,
    50: 57,
    51: 58,
    52: 59,
    53: 60,
    54: 61,
    55: 63,
    56: 64,
    57: 65,
    58: 66,
    59: 67,
    60: 68,
    61: 69,
    62: 72,
    63: 73,
    64: 74,
    65: 75,
    66: 76,
    67: 77,
    68: 78,
    69: 79,
    70: 80,
    71: 81,
    72: 82,
    73: 84,
    74: 85,
    75: 86,
    76: 87,
    77: 88,
    78: 89,
    79: 90,
    80: 91,
    81: 92,
    82: 93,
    83: 94,
    84: 95,
    85: 96,
    86: 97,
    87: 98,
    88: 99,
    89: 100,
    90: 101,
    91: 102,
    92: 103,
}

_index_to_imgt_state[("imgt", "GA1")] = _index_to_imgt_state[("imgt", "MH1")]
_index_to_imgt_state[("imgt", "GA2")] = dict(
    [(k + 1000, v) for k, v in list(_index_to_imgt_state[("imgt", "GA1")].items())]
)


# To map indices onto index_to_imgt_state: this works for everything apart from MHC constant domains.
def map_state_index_imgt(statevector):
    numdict = {}
    curr_num = 1
    for i, char in enumerate(statevector):
        if not char.isdigit():
            numdict[curr_num] = i
            curr_num += 1
    return numdict


# To map indices onto index_to_imgt_state
def map_state_index_mhc_clike(statevector):
    numdict = {}
    curr_num = 1
    for i, char in enumerate(statevector):
        if not char.isdigit():
            numdict[curr_num] = i
            if curr_num == 45:
                curr_num += 32
            else:
                curr_num += 1
    return numdict


_index_to_imgt_state[("imgt", "B2M")] = map_state_index_mhc_clike(mhc_cvec)

# For N numbering, IMGT uses numbers 1000-1100.
# For now, we'll avoid numbering B2M domains and the C-LIKE domains.
_index_to_imgt_state[("imgt", "MH1")].update(
    dict(
        [(k + 1000, v) for k, v in list(_index_to_imgt_state[("imgt", "MH1")].items())]
    )
)
_index_to_imgt_state[("imgt", "GA")] = _index_to_imgt_state[("imgt", "GA1")]
_index_to_imgt_state[("imgt", "GB")] = _index_to_imgt_state[("imgt", "GA1")]
_index_to_imgt_state[("imgt", "GA1L")] = _index_to_imgt_state[("imgt", "CD1")]
_index_to_imgt_state[("imgt", "GA2L")] = dict(
    [(k + 1000, v) for k, v in list(_index_to_imgt_state[("imgt", "CD1")].items())]
)
_index_to_imgt_state[("imgt", "CD1")].update(
    dict(
        [(k + 1000, v) for k, v in list(_index_to_imgt_state[("imgt", "CD1")].items())]
    )
)
_index_to_imgt_state[("imgt", "MR1")] = _index_to_imgt_state[("imgt", "CD1")]


def get_region(position, chain):
    """
    Get the region in which the position belongs given the chain, numbering scheme and definition.

    **Note** this function does not know about insertions on the sequence. Therefore, it will get the region annotation
    wrong when using non-equivalent scheme-definitions.

    To get around this please use the annotate_regions function which implements heuristics to get the definition correct
    in the scheme.

    """
    index, insertion = position
    chain = chain.upper()

    # imgt_state is a dictionary that maps an IMGT position (e.g. 1) onto the position along the state vector;
    # Thus, 1 maps to 0 because 1 is the first number in the IMGT numbering scheme but is the 0th position along the state vector.
    imgt_state = _index_to_imgt_state[("imgt", chain)]
    # Get the state vector corresponding to a particular chain; this is either 1,2,3,4,5,6,7, corresponding to different TCR regions.
    state_vec = _regions["imgt"][chain]

    if chain in TCR_CHAINS:
        if index in imgt_state:
            state_idx = imgt_state[index]
            # Returns a fwb1, cdra2... etc.
            return _reg_one2three[state_vec[state_idx]] % chain.lower()
        else:
            return "?"
    else:
        if index in imgt_state:
            state_idx = imgt_state[index]
            # Returns whether helix or turn on the MHC (G-Domain annotation)
            return _reg_tostring[state_vec[state_idx]]
        else:
            return "?"


def annotate_regions(numbered_sequence, chain):
    """
    Given a numbered sequence (list) annotate which region each residue belongs to.
    Currently, only the IMGT numbering and definition are implemented.
    If possible, use the corresponding numbering scheme and definition.

    This function automates the heuristics recognise different definitions in each scheme. However,
    some of the conversions are non-trivial.
    """
    chain = chain.upper()
    c = chain.lower()

    numdict = dict(numbered_sequence)

    cdr_acceptors = {1: Accept(), 2: Accept(), 3: Accept()}

    cdr_acceptors[1].set_regions(["cdr%s1" % c])
    cdr_acceptors[2].set_regions(["cdr%s2" % c])
    cdr_acceptors[3].set_regions(["cdr%s3" % c])

    # We start off by annotating framework regions; this switches when we find a CDR
    fw_region = "fw%s1" % c
    region_annotations = []

    cterm = max(_index_to_imgt_state[("imgt", chain)].keys())
    for r, a in numbered_sequence:
        if cdr_acceptors[1].accept(r, chain):
            region_annotations.append((r, a, "cdr%s1" % c))
            fw_region = "fw%s2" % c
        elif cdr_acceptors[2].accept(r, chain):
            region_annotations.append((r, a, "cdr%s2" % c))
            fw_region = "fw%s3" % c
        elif cdr_acceptors[3].accept(r, chain):
            region_annotations.append((r, a, "cdr%s3" % c))
            fw_region = "fw%s4" % c
        elif (
            r[0] <= cterm
        ):  # Anything out of the variable region is not assigned a region i.e. ''
            region_annotations.append((r, a, fw_region))
        else:
            region_annotations.append((r, a, ""))

    return region_annotations


class Accept(object):
    """
    A class to select which positions should be compared.
    """

    _defined_regions = TCR_REGIONS
    _macro_regions = {
        "bframework": set(["fwb1", "fwb2", "fwb3", "fwb4"]),
        "bcdrs": set(["cdrb1", "cdrb2", "cdrb3"]),
        "aframework": set(["fwa1", "fwa2", "fwa3", "fwa4"]),
        "acdrs": set(["cdra1", "cdra2", "cdra3"]),
        "gframework": set(["fwg1", "fwg2", "fwg3", "fwg4"]),
        "gcdrs": set(["cdrg1", "cdrg2", "cdrg3"]),
        "dframework": set(["fwd1", "fwd2", "fwd3", "fwd4"]),
        "dcdrs": set(["cdrd1", "cdrd2", "cdrd3"]),
    }

    _macro_regions.update(
        {
            "framework": _macro_regions["bframework"]
            | _macro_regions["aframework"]
            | _macro_regions["dframework"]
            | _macro_regions["gframework"],
            "cdrs": _macro_regions["bcdrs"]
            | _macro_regions["acdrs"]
            | _macro_regions["dcdrs"]
            | _macro_regions["gcdrs"],
            "vb": _macro_regions["bcdrs"] | _macro_regions["bframework"],
            "va": _macro_regions["acdrs"] | _macro_regions["aframework"],
            "vg": _macro_regions["gcdrs"] | _macro_regions["gframework"],
            "vd": _macro_regions["dcdrs"] | _macro_regions["dframework"],
        }
    )

    _macro_regions.update(
        {
            "ba": _macro_regions["vb"] | _macro_regions["va"],
            "dg": _macro_regions["vd"] | _macro_regions["vg"],
            "tr": _macro_regions["vb"]
            | _macro_regions["va"]
            | _macro_regions["vd"]
            | _macro_regions["vg"],
        }
    )

    _macro_positions = {}

    def __init__(self, NOT=False):
        self.NOT = NOT
        self.set_regions()
        self.positions = {"B": set(), "A": set(), "D": set(), "G": set()}
        self.exclude = {"B": set(), "A": set(), "D": set(), "G": set()}

    def set_regions(self, regions=[]):
        """
        Set the regions to be used. Will clear anything added using add regions.
        """
        if self.NOT:
            self.regions = self._macro_regions["tr"]
        else:
            self.regions = set()
        self.add_regions(regions)

    def add_regions(self, regions):
        """
        Add regions to the selection.
        """
        for region in regions:
            region = region.lower()
            if region in self._defined_regions:
                if self.NOT:
                    self.regions = self.regions - set([region])
                else:
                    self.regions.add(region)
            elif region in self._macro_regions:
                if self.NOT:
                    self.regions = self.regions - self._macro_regions[region]
                else:
                    self.regions = self.regions | self._macro_regions[region]
            elif region in self._macro_positions:
                raise AssertionError("Undefined region")
            else:
                raise AssertionError("Undefined region")

    def add_positions(self, positions, chain):
        for position in positions:
            index, insertion = position
            self.positions[chain].add((index, insertion))

    def exclude_positions(self, positions, chain):
        for position in positions:
            index, insertion = position
            self.exclude[chain].add((index, insertion))

    def accept(self, position, chain):
        if position in self.exclude[chain]:
            return 0
        elif (
            get_region(position, chain) in self.regions
            or position in self.positions[chain]
        ):
            return 1
        return 0
