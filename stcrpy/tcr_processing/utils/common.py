import math
import warnings
import numpy as np


def fastcross(v, w):
    """Cross-vector of two Vector objects which is faster than NumPy's version"""
    return np.array(
        [
            v[1] * w[2] - v[2] * w[1],
            v[2] * w[0] - v[0] * w[2],
            v[0] * w[1] - v[1] * w[0],
        ]
    )


def fastnorm(A):
    """Faster version of Euclidean norm"""
    return math.sqrt(sum([x**2 for x in A]))


def identity(seq1, seq2, positions=[]):
    """
    Find the matched sequence identity between two aligned sequences.
    Can accept lists/strings, but this assumes that the two sequences are of the same length.
    Args:
        seq1: Dictionary with key as the position and value as the single letter amino acid code. or an aligned list or string
        seq2: Dictionary with key as the position and value as the single letter amino acid code. or an aligned list or string
    """
    n = 0  # number
    m = 0  # match

    if isinstance(seq1, dict) and isinstance(seq2, dict):
        if not positions:
            positions = set(seq1.keys()) | set(seq2.keys())
    else:
        assert len(seq1) == len(seq2), "Use two aligned sequences."
        positions = range(len(seq1))

    # matched identity
    for p in positions:
        try:
            if seq1[p] == "-":
                continue
            if seq2[p] == "-":
                continue
        except KeyError:
            continue

        if seq1[p] == seq2[p]:
            m += 1
        n += 1

    try:
        return float(m) / n
    except ZeroDivisionError:
        return 0


def angle(v1, v2):
    """Return the angle between two vectors"""
    # num = np.dot(v1.point,v2.point)
    # denom = v1.norm() * v2.norm()
    num = np.dot(v1, v2)
    denom = fastnorm(v1) * fastnorm(v2)
    if abs(num / denom) > 1:
        return np.pi
    return np.arccos(num / denom)
