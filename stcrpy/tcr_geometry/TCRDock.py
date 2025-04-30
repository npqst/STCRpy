#!/usr/bin/env python

"""
__init__.py
Description: 	Calculate the TCR/pMHC docking angle.
                Based on the method by Rudolph, Stanfield, Wilson (Annu Revl Immmunol 2006, 24:419-466).
                This method is preferred for cases when the antigen is not a peptide.
                Formula uses SVD.
Jun 12, 2017
"""
# from TCRDB.TcrPDB.TCR import TCR
from ..tcr_processing import TCR, MHC, MHCchain
from ..utils.error_stream import ErrorStream

import sys
import warnings
import numpy as np


class TCRDock(object):
    def __init__(self, tcr, QUIET=False):
        """
        Calculate the docking angle between TCR and pMHC.

        Args:
            tcr (TCR): TCR input
            QUIET (bool, optional): Verbosity of error stream. Defaults to False.
        """
        self.warnings = ErrorStream()
        self.QUIET = QUIET

        self.TCR = tcr
        self.angle = np.nan

        # Get the MHC for the TCR.
        self.MHC = tcr.get_MHC()

        if not self.MHC:
            self.warnings.write(
                "The TCR structure does not have a detected MHC molecule. No docking angle will be calculated.\n"
            )
            self.abort = True
            return

        self.MHC = self.MHC[0]

        self.abort = False
        if not isinstance(self.TCR, TCR):
            self.warnings.write(
                "The TCR structure is an unpaired TCR chain. No docking angle will be calculated.\n"
            )
            self.abort = True
            return

        if isinstance(self.MHC, MHC):
            pass

        elif not isinstance(self.MHC, MHC) and isinstance(self.MHC, MHCchain):
            if (
                self.MHC.chain_type == "MH1"
                or self.MHC.chain_type == "CD1"
                or self.MHC.chain_type == "MR1"
                or self.MHC.chain_type == "GA1"
            ):
                acceptable_range = list(range(50, 87))
                residues = [
                    r
                    for r in self.MHC.get_residues()
                    if r.id[1] % 1000 in acceptable_range
                ]
                if len(residues) >= (len(acceptable_range) - 10):
                    self.warnings.write(
                        "Warning: detected an MHC chain of type %s; doesn't seem to have an associated B2M molecule.\n"
                        % self.MHC.chain_type
                    )
                    pass
                else:
                    self.warnings.write(
                        "An MHC molecule was not found. No docking angle will be calculated.\n"
                    )
                    self.abort = True
                    return
            else:
                self.warnings.write(
                    "An MHC molecule was not found. No docking angle will be calculated.\n"
                )
                self.abort = True
                return

        elif not isinstance(self.MHC, MHC):
            self.warnings.write(
                "An MHC molecule was not found. No docking angle will be calculated.\n"
            )
            self.abort = True
            return

        # Resolve the vectors for the TCR
        self._resolve_vectors()

    def _resolve_vectors(self):
        if self.abort:
            return
        # Get the vector between cysteine centroids.
        self._get_cysteine_vector()

        # Get the vector of helices
        self._get_helix_vectors()

    def _get_cysteine_vector(self):
        """
        Get the centroids of the disulphide bridge and calculate a vector through it.
        """
        # Get variable domains
        if self.TCR.get_TCR_type() == "abTCR":
            vbg, vda = self.TCR.get_VB(), self.TCR.get_VA()
        elif self.TCR.get_TCR_type() == "gdTCR":
            vbg, vda = self.TCR.get_VD(), self.TCR.get_VG()
        elif self.TCR.get_TCR_type() == "dbTCR":
            vbg, vda = self.TCR.get_VB(), self.TCR.get_VD()

        try:
            # Get sulphur atoms of each of the cysteines
            bg_23, bg_104 = vbg[23]["SG"], vbg[104]["SG"]
            da_23, da_104 = vda[23]["SG"], vda[104]["SG"]
            bg_centroid = np.mean((bg_23.coord, bg_104.coord), axis=0)
            da_centroid = np.mean((da_23.coord, da_104.coord), axis=0)

            # Compute the vector between the centroids
            self.vec_centroid = bg_centroid - da_centroid

        except KeyError:
            self.warnings.write(
                "Cysteine(s) or sulphur atom(s) not detected. Check for IMGT residues 23/104 in beta/alpha/delta/gamma chains.\n"
            )
            self.abort = True
            return

    def _get_helix_vectors(self):
        """
        Get the best fit vector for the CA atoms.
        For MH1 and MH2, the atoms are based on the positions from Rudolph et al., 2006 with the IMGT numbering;
        For CD1 and MR1, we use the same rules as MH1.
        """
        try:
            if self.MHC.get_MHC_type() == "MH1":

                # Get CA atoms of 50-86 and 1050-1086 (A140-A176 on Rudolph et al).
                # Using the modulus operator helps to get the last 2 digits of the IMGT-numbered residue. https://stackoverflow.com/a/28570538
                acceptable_range = list(range(50, 87))
                ca_atoms = np.array(
                    [
                        r["CA"].coord
                        for r in self.MHC.get_alpha().get_residues()
                        if r.id[1] % 1000 in acceptable_range and "CA" in r
                    ]
                )

            elif self.MHC.get_MHC_type() == "CD1":
                # Get CA atoms of 50-86 and 1050-1086 (A140-A176 on Rudolph et al).
                # Using the modulus operator helps to get the last 2 digits of the IMGT-numbered residue. https://stackoverflow.com/a/28570538
                acceptable_range = list(range(50, 87))
                ca_atoms = np.array(
                    [
                        r["CA"].coord
                        for r in self.MHC.get_CD1().get_residues()
                        if r.id[1] % 1000 in acceptable_range and "CA" in r
                    ]
                )

            elif self.MHC.get_MHC_type() == "MR1":
                # Get CA atoms of 50-86 and 1050-1086 (A140-A176 on Rudolph et al).
                # Using the modulus operator helps to get the last 2 digits of the IMGT-numbered residue. https://stackoverflow.com/a/28570538
                acceptable_range = list(range(50, 87))
                ca_atoms = np.array(
                    [
                        r["CA"].coord
                        for r in self.MHC.get_MR1().get_residues()
                        if r.id[1] % 1000 in acceptable_range and "CA" in r
                    ]
                )

            elif self.MHC.get_MHC_type() == "MH2":
                # Get CA atoms of A and B52-87
                # Using the modulus operator helps to get the last 2 digits of the IMGT-numbered residue. https://stackoverflow.com/a/28570538
                alpha_range = list(range(50, 88))
                beta_range = alpha_range[2:]
                ca_atoms = [
                    r["CA"].coord
                    for r in self.MHC.get_GA().get_residues()
                    if r.id[1] in alpha_range and "CA" in r
                ]
                ca_atoms += [
                    r["CA"].coord
                    for r in self.MHC.get_GB().get_residues()
                    if r.id[1] in beta_range and "CA" in r
                ]
                ca_atoms = np.array(ca_atoms)

        except AttributeError:
            if (
                self.MHC.chain_type == "MH1"
                or self.MHC.chain_type == "CD1"
                or self.MHC.chain_type == "MR1"
                or self.MHC.chain_type == "GA1"
            ):
                acceptable_range = list(range(50, 87))
                ca_atoms = np.array(
                    [
                        r["CA"].coord
                        for r in self.MHC.get_residues()
                        if r.id[1] % 1000 in acceptable_range and "CA" in r
                    ]
                )
            else:
                self.abort = True
                return

        self.ca = ca_atoms

    def calculate_docking_angle(self, force=False):
        if not np.isnan(self.angle):
            return self.angle
        elif force:
            self._resolve_vectors()
        elif self.abort:
            return np.nan

        # Compute the mean and calculate the vector using SVD
        # https://stackoverflow.com/a/2333251
        ca_centroid = self.ca.mean(axis=0)
        centred_dat = self.ca - ca_centroid
        u, d, v = np.linalg.svd(centred_dat)

        # The first row of v is the 1st principal component.
        self.V = v[0]

        self.angle = self._angle(self.V, self.vec_centroid)

        if not self.QUIET and self.warnings.log:
            sys.stderr.write("\n".join(self.warnings.log))
            sys.stderr.write("\n")

        return self.angle

    def _angle(self, v1, v2):
        """
        Return the angle between two vectors in degrees.
        print an error message if the numerator is negative
        """
        # Check the direction of the dot product; assert positive, as we know the angle should be between 0-90 deg.
        # This is because the singular-value decomposition for finding the best fit might return a different sign than we require
        # https://stackoverflow.com/questions/17682626/singular-value-decomposition-different-results-with-jama-pcolt-and-numpy
        # https://math.stackexchange.com/questions/2359992/how-to-resolve-the-sign-issue-in-a-svd-problem

        numerator = np.dot(v1, v2)
        denominator = np.linalg.norm(v1) * np.linalg.norm(v2)

        if numerator < 0:
            numerator = abs(numerator)

        if numerator / denominator > 1.0:
            return 180.0
        else:
            return np.degrees(np.arccos(numerator / denominator))
