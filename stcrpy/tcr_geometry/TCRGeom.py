import warnings
import numpy as np

from ..tcr_processing import abTCR, MHCchain
from .TCRCoM import MHCI_TCRCoM, MHCII_TCRCoM


class TCRGeom:
    """Class for TCR geometry calculations."""

    def __init__(
        self,
        tcr: "TCR",
        save_aligned_as: bool = False,
        polarity_as_sign: bool = True,
        mode: str = "cys",
    ):
        """Constructor for TCR geometry calculation class.

        Args:
            tcr (TCR): TCR structure object of which to calculate complex geometry. TCR must be in complex to MHC.
            save_aligned_as (bool, optional): Save the alignment of the TCR:pMHC complex as PDB file for verification. Defaults to False.
            polarity_as_sign (bool, optional): Use the sign of the angle to indicate non-canonical pose. ie reverse binding angles are assigned negative sign. Defaults to True.
            mode (str, optional): Method to calculate geometry with, see STCRpy paper for details. Defaults to "cys". Options: 'rudolph', 'com', 'cys'.
        """
        self.tcr_com = None
        self.mhc_com = None
        self.tcr_VA_com = None
        self.tcr_VB_com = None
        self.tcr_VA_cys_centroid = None
        self.tcr_VB_cys_centroid = None
        self.tcr_VA_VB_angle = None
        self.scanning_angle = np.nan
        self.tcr_pitch_angle = np.nan
        self.tcr_docking_angles_cys = None
        self.tcr_mhc_dist = None
        self.mode = mode

        self._type_checks(tcr)

        mhc = tcr.get_MHC()[0]

        self._set_mhc_reference(mhc)

        if self.mode != "rudolph":
            (self.tcr_com, self.mhc_com, self.tcr_VA_com, self.tcr_VB_com) = (
                self.mhc_tcr_com_calculator.calculate_centres_of_mass(
                    tcr, save_aligned_as=save_aligned_as
                )
            )

        if self.mode != "com":
            (self.tcr_VA_cys_centroid, self.tcr_VB_cys_centroid) = (
                self._get_cys_centroids(tcr)
            )

        if self.mode in ["cys", "rudolph"]:
            self.tcr_vector = self.get_tcr_vector(
                self.tcr_VA_cys_centroid, self.tcr_VB_cys_centroid
            )
        elif self.mode == "com":
            self.tcr_vector = self.get_tcr_vector(self.tcr_VA_com, self.tcr_VB_com)
        else:
            warnings.warn(
                f"Geometry mode {self.mode} not recognised. Using CYS atom coordinates to calculate geometry"
            )
            self.tcr_vector = self.get_tcr_vector(
                self.tcr_VA_cys_centroid, self.tcr_VB_cys_centroid
            )

        if self.mode == "rudolph":
            self.mhc_vector = self._get_mhc_helix_vectors(mhc)
            self.scanning_angle = self.calculate_rudolph_angle(
                self.tcr_vector, self.mhc_vector
            )

        else:
            self.scanning_angle, self.tcr_pitch_angle = (
                self.calculate_tcr_docking_angles(
                    self.tcr_vector, polarity_as_sign=polarity_as_sign
                )
            )

        self.polarity = self.get_polarity()

    def __repr__(self) -> str:
        """Print TCR:pMHC complex geometry.

        Returns:
            str: TCR:pMHC complex geometry metrics.
        """
        def polarity_to_str(polarity):
            return "canonical" if polarity == 0 else "reverse"

        def mode_to_str(mode):
            if mode == "cys":
                return "Cysteine centroids"
            elif mode == "com":
                return "Centre of mass"
            elif mode == "rudolph":
                return "Rudolph et. al. 2006"
            else:
                return mode

        return (
            f"TCR CoM: {self.tcr_com}\nMHC CoM: {self.mhc_com}\n"
            + f"TCR VA CoM: {self.tcr_VA_com}\nTCR VB CoM: {self.tcr_VB_com}\n"
            + f"TCR VA CYS Centroid: {self.tcr_VA_cys_centroid}\nTCR VB CYS Centroid: {self.tcr_VB_cys_centroid}\n"
            + f"Scanning angle: {np.degrees(self.scanning_angle)}\n"
            + f"Pitch angle: {np.degrees(self.tcr_pitch_angle)}\n"
            + f"Polarity: {polarity_to_str(self.polarity)}\n"
            + f"Geometry mode: {mode_to_str(self.mode)}"
        )

    def to_dict(self) -> dict:
        """Return TCR:pMHC complex geometry metrics as dictionary.

        Returns:
            dict: TCR:pMHC complex geometry.
        """
        return {
            "tcr_com": [self.tcr_com],
            "mhc_com": [self.mhc_com],
            "tcr_VA_com": [self.tcr_VA_com],
            "tcr_VB_com": [self.tcr_VB_com],
            "tcr_VA_cys_centroid": [self.tcr_VA_cys_centroid],
            "tcr_VB_cys_centroid": [self.tcr_VB_cys_centroid],
            "scanning_angle": np.degrees(self.scanning_angle),
            "pitch_angle": np.degrees(self.tcr_pitch_angle),
            "polarity": self.polarity,
            "mode": self.mode,
        }

    def to_df(self) -> "pandas.DataFrame":
        """Return TCR:pMHC complex geometry metrics and pandas dataframe.

        Returns:
            pandas.DataFrame: TCR:pMHC complex geometry metrics
        """
        import pandas as pd

        return pd.DataFrame.from_dict(self.to_dict())

    def get_scanning_angle(self, rad: bool = False) -> float:
        """Return TCR:pMHC complex scanning (aka crossing, incident angle) of TCR to MHC.

        Args:
            rad (bool, optional): Return angle in radians. Defaults to False.

        Returns:
            float: TCR:pMHC scanning angle.
        """
        if rad:
            return self.scanning_angle
        else:
            return np.degrees(self.scanning_angle)

    def get_pitch_angle(self, rad: bool = False) -> float:
        """Return TCR:pMHC pitch angle, ie tilt of the TCR over the MHC.

        Args:
            rad (bool, optional): Return angle in radians. Defaults to False.

        Returns:
            float: TCR:pMHC pitch angle
        """
        if rad:
            return self.tcr_pitch_angle
        else:
            return np.degrees(self.tcr_pitch_angle)

    def _type_checks(self, tcr: "TCR"):
        """Run checks on TCR structure object to ensure geometry can be calculated.
        \nChecks if TCR is in complex with MHC and whether TCR is alpha/beta chain TCR.

        Args:
            tcr (TCR): TCR structure object.

        Raises:
            NotImplementedError: TCR:pMHC geometry calculations are compatible the alpha beta chain TCRs.
        """
        if len(tcr.get_MHC()) == 0:
            warnings.warn(
                f"No MHC associated with TCR {tcr.parent.parent.id}_{tcr.id}. Geometry cannot be calculated."
            )
            return

        if not isinstance(tcr, abTCR):
            raise NotImplementedError(
                f"TCR MHC geometry only implemented for abTCR types, not {type(tcr)}"
            )

    def _set_mhc_reference(self, mhc: "MHC"):
        """Sets the MHC structure reference to use for alignment.

        Args:
            mhc (MHC): MHC structure object

        Raises:
            NotImplementedError: MHC alignments are compatible with MHC class I and II. CD1 and MR1 types currently not supported.
            ValueError: MHC unrecognised.
        """
        if (isinstance(mhc, MHCchain) and mhc.chain_type not in ["GA", "GB"]) or (
            hasattr(mhc, "MHC_type") and mhc.MHC_type == "MH1"
        ):
            self.mhc_tcr_com_calculator = MHCI_TCRCoM()
        elif (isinstance(mhc, MHCchain) and mhc.chain_type in ["GA", "GB"]) or (
            hasattr(mhc, "MHC_type") and mhc.MHC_type == "MH2"
        ):
            self.mhc_tcr_com_calculator = MHCII_TCRCoM()
        else:
            if hasattr(mhc, "MHC_type") and mhc.get_MHC_type() in ["CD1", "MR1"]:
                raise NotImplementedError(
                    "TCR geometry not yet implemented for CD1 and MR1 antigen."
                )
            else:
                raise ValueError(f"MHC type of {mhc} not recognised.")

    def get_tcr_vector(self, tcr_VA_com: np.array, tcr_VB_com: np.array) -> np.array:
        """Calculates the vector from the VA centre of mass to the VB centre of mass, translated such that the vector originates at the total TCR's centre of mass and truncated to unit length.
        Args:
            tcr_VA_com (np.array): TCR VA centre of mass
            tcr_VB_com (np.array): TCR VB center of mass

        Returns:
            np.array: Unit vector from TCR CoM in direction of VA CoM to VB CoM
        """
        direction_vec = tcr_VB_com - tcr_VA_com
        tcr_vector = direction_vec / np.linalg.norm(direction_vec)
        return tcr_vector

    def calculate_tcr_docking_angles(
        self, tcr_vector: np.array, polarity_as_sign: bool = True
    ) -> tuple[np.array]:
        """
        Calculates the scanning angle and pitch of the TCR relative to the MHC.
        This function relies on the previous alignment of the MHC to the reference MHC,
        with the x axis defined perpedicular to the peptide, the y-axis along the peptide,
        and the z-axis pointing 'up' away from the MHC.

        The scanning angle is calculated as the angle between the y axis
        and the projection of the TCR vector, which points from VA to VB, onto the x-y plane.

        The pitch is calculated as the angle between the z-axis and the plane normal
        to the TCR vector and passing through the TCR CoM. Specifically the angle of the vector from the
        z-axis intersection and the nearest projection of the z-axis onto the plane dividing the VA and VB
        TCR domains is calculated, which is defined as the pitch of the TCR.

        Args:
            tcr_vector (np.array): Unit vector in direction of VA CoM to VB CoM
            polarity (bool, optional): Set the polarity as the sign of the scanning angle
                                        ie negative if polarity is reversed (1),
                                        else positive if polarity is canonical (0). Defaults to True.
        Returns:
            tuple[np.array]: Tuple containing the scanning angle and pitch angle of the TCR to the MHC
        """

        xy_projection = tcr_vector[:2] / np.linalg.norm(tcr_vector[:2])
        self.scanning_angle = np.arccos(np.dot(xy_projection, np.asarray([0.0, 1.0])))
        phi = np.arccos(np.sqrt(1 - (tcr_vector[-1] ** 2)))
        if polarity_as_sign:
            self.scanning_angle = self.scanning_angle * ((-1) ** self.get_polarity())
        return self.scanning_angle, phi

    def get_polarity(self) -> int:
        """
        Return the polarity of the TCR based on the TCR vector pointing from the VA to the VB CoM, or the scanning angle if using Rudolph et. al..
        If the x component is negative, ie the tcr_vector points from the alpha 2 chain to the
        alpha 1 chain of the MHC, the polarity is canonical (0). Otherwise the polarity is reverse (1).

        Returns:
            int: 0 for canonical polarity, 1 for reverse polarity.
        """
        if self.mode in ["cys", "com"]:
            if self.tcr_vector[0] <= 0:
                return 0
            else:
                return 1
        elif self.mode == "rudolph":
            return 0 if np.abs(np.degrees(self.scanning_angle)) < 120.0 else 1

    def _get_cys_centroids(self, tcr: "TCR") -> tuple[np.array, np.array]:
        """Calculate the halfway coordintate bewtween the CYS residues of the V.

        Args:
            tcr (TCR): TCR structure object

        Raises:
            KeyError: Cys residue not found in V domain
            KeyError: SG and CA atoms not found in Cys residue

        Returns:
            tuple[np.array, np.array]: (Cys centroid of VA/VG, Cys centroid of VB/VD)
        """
        domain_order = {
            "VA": 0,
            "VB": 3,
            "VG": 1,
            "VD": 2,
        }  # defines order s.t. order is VA->VB, VG->VD, VD->VB depending on pairing
        cys_coords = {dom: {} for dom in tcr.get_domain_assignment()}
        for domain, chain in tcr.get_domain_assignment().items():
            for cys_res_nr in [23, 104]:
                try:
                    cys_coords[domain][cys_res_nr] = tcr[chain][cys_res_nr]["SG"].coord
                except KeyError as e:
                    if str(cys_res_nr) in str(e):
                        raise KeyError(
                            f"IMGT numbered residue {str(cys_res_nr)} not found in {tcr.id} domain {domain} with chain ID {chain}. Consider calculating TCR geometry with centre of mass coordinates. "
                        )
                    elif "SG" in str(e):
                        warnings.warn(
                            f"SG atom not found in IMGT residue number {str(cys_res_nr)} of {tcr.id} domain {domain} with chain ID {chain}. Trying CA atom instead."
                        )
                        try:
                            cys_coords[domain][cys_res_nr] = tcr[chain][cys_res_nr][
                                "CA"
                            ]
                        except KeyError:
                            raise KeyError(
                                f"Neither SG not CA atom found in IMGT residue number {str(cys_res_nr)} of {tcr.id} domain {domain} with chain ID {chain}. Consider calculating TCR geometry with centre of mass coordinates."
                            )

        cys_centroids = (
            np.mean([cys_coords[dom][res_nr] for res_nr in [23, 104]], axis=0)
            for dom in sorted(
                cys_coords, key=lambda x: domain_order[x]
            )  # calculate cys centroids and orders such that vector can be calculated from VA to VB (VG to VD)
        )

        return cys_centroids  # VA CYS centroid, VB CYS centroid

    def _get_mhc_helix_vectors(self, mhc: "MHC") -> np.array:
        """Calculate the vector pointing along the MHC using the MHC helices. Points approximately from N to C terminus of GA chain.
        \nVector is calculated as principal component of SVD decomposition of point cloud of backbone carbon alpha atoms in the helices forming the peptide cleft.

        Args:
            mhc (MHC): MHC structure object

        Returns:
            np.array: Unit vector along the MHC
        """
        helix_residue_ranges = {
            "MH1": ((50, 87), (140, 177)),
            "MH2": {"GA": ((50, 88),), "GB": ((54, 65), (67, 91))},
            # "CD1": (50, 87),              # TODO: determine consisent residue ranges for CD1 and MR1 helices
            # "MR1": (50, 87),
        }

        if mhc.get_MHC_type() == "MH1":
            helix_CA_coords = np.asarray(
                [
                    mhc.get_MH1()[(" ", i, " ")]["CA"].coord
                    for _range in helix_residue_ranges[mhc.get_MHC_type()]
                    for i in range(*_range)
                    if (" ", i, " ") in mhc.get_MH1()
                    and "CA" in mhc.get_MH1()[(" ", i, " ")]
                ]
            )
            # mhc vector should point approximately from N to C terminii of GA1 (MH1)
            try:  # wrap in try except in case residues at edges don't exist or atoms are missing
                approximate_mhc_vector = (
                    mhc.get_MH1()[(" ", helix_residue_ranges["MH1"][0][1] - 1, " ")][
                        "CA"
                    ].coord
                    - mhc.get_MH1()[(" ", helix_residue_ranges["MH1"][0][0], " ")][
                        "CA"
                    ].coord
                )
            except KeyError:  # try surrounding residues
                approximate_mhc_vector = (
                    mhc.get_MH1()[(" ", helix_residue_ranges["MH1"][0][1] - 2, " ")][
                        "CA"
                    ].coord
                    - mhc.get_MH1()[(" ", helix_residue_ranges["MH1"][0][0] + 1, " ")][
                        "CA"
                    ].coord
                )
            approximate_mhc_vector = approximate_mhc_vector / np.linalg.norm(
                approximate_mhc_vector
            )
        elif mhc.get_MHC_type() == "MH2":
            helix_CA_coords = np.asarray(
                [
                    mhc.get_GA()[(" ", i, " ")]["CA"].coord
                    for _range in helix_residue_ranges[mhc.get_MHC_type()]["GA"]
                    for i in range(*_range)
                    if (" ", i, " ") in mhc.get_GA()
                    and "CA" in mhc.get_GA()[(" ", i, " ")]
                ]
                + [
                    mhc.get_GB()[(" ", i, " ")]["CA"].coord
                    for _range in helix_residue_ranges[mhc.get_MHC_type()]["GB"]
                    for i in range(*_range)
                    if (" ", i, " ") in mhc.get_GB()
                    and "CA" in mhc.get_GB()[(" ", i, " ")]
                ]
            )
            # mhc vector should point approximately from N to C terminii of GA (MH2)
            try:  # wrap in try except in case residues at edges don't exist or atoms are missing
                approximate_mhc_vector = (
                    mhc.get_GA()[
                        (" ", helix_residue_ranges["MH2"]["GA"][0][1] - 1, " ")
                    ]["CA"].coord
                    - mhc.get_GA()[(" ", helix_residue_ranges["MH2"]["GA"][0][0], " ")][
                        "CA"
                    ].coord
                )
            except KeyError:  # try surrounding residues
                approximate_mhc_vector = (
                    mhc.get_GB()[
                        (" ", helix_residue_ranges["MH2"]["GA"][0][1] - 2, " ")
                    ]["CA"].coord
                    - mhc.get_GB()[
                        (" ", helix_residue_ranges["MH2"]["GA"][0][0] + 1, " ")
                    ]["CA"].coord
                )
            approximate_mhc_vector = approximate_mhc_vector / np.linalg.norm(
                approximate_mhc_vector
            )

        mhc_helix_points = helix_CA_coords - np.mean(
            helix_CA_coords, axis=0
        )  # centre helices
        mhc_vector = np.linalg.svd(mhc_helix_points)[2][0]  # fit line
        # check if approximate vector and MHC vector are aligned by checking magnitude of vector sum, if not, invert MHC vector
        if np.linalg.norm(mhc_vector + approximate_mhc_vector) < np.linalg.norm(
            mhc_vector - approximate_mhc_vector
        ):
            mhc_vector = -mhc_vector
        mhc_vector = mhc_vector / np.linalg.norm(mhc_vector)
        return mhc_vector

    def calculate_rudolph_angle(
        self, tcr_vector: np.array, mhc_vector: np.array
    ) -> float:
        """Calculate the scanning angle of TCR:pMHC complex according to Rudolph et. al.

        Args:
            tcr_vector (np.array): Vector pointing from VA Cys centroid to VB Cys centroid
            mhc_vector (np.array): Vector pointing along MHC cleft

        Returns:
            float: scanning angle
        """
        scanning_angle = np.arccos(
            np.dot(tcr_vector, mhc_vector)
            / (np.linalg.norm(tcr_vector) * np.linalg.norm(mhc_vector))
        )
        return scanning_angle
