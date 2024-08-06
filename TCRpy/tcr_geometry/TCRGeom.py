import warnings
import numpy as np

from ..tcr_processing import abTCR, MHCchain
from .TCRCoM import MHCI_TCRCoM, MHCII_TCRCoM


class TCRGeom:
    def __init__(self, tcr, save_aligned_as=False, polarity_as_sign=True):
        self.tcr_com = None
        self.mhc_com = None
        self.tcr_VA_com = None
        self.tcr_VB_com = None
        self.tcr_VA_cys_coords = None
        self.tcr_VB_cys_coords = None
        self.tcr_VA_VB_angle = None
        self.tcr_scanning_angle = None
        self.tcr_pitch_angle = None
        self.tcr_docking_angles_cys = None
        self.tcr_mhc_dist = None

        if len(tcr.get_MHC()) == 0:
            warnings.warn("No MHC associated with TCR. Geometry cannot be calculated.")
        mhc = tcr.get_MHC()
        if not isinstance(tcr, abTCR):
            raise NotImplementedError(
                f"TCR MHC geometry only implemented for abTCR types, not {type(tcr)}"
            )

        if (isinstance(mhc[0], MHCchain) and mhc[0].chain_type not in ["GA", "GB"]) or (
            hasattr(mhc[0], "MHC_type") and mhc[0].MHC_type == "MH1"
        ):
            self.mhc_tcr_com_calculator = MHCI_TCRCoM()
        elif (isinstance(mhc[0], MHCchain) and mhc[0].chain_type in ["GA", "GB"]) or (
            hasattr(mhc[0], "MHC_type") and mhc[0].MHC_type == "MH2"
        ):
            self.mhc_tcr_com_calculator = MHCII_TCRCoM()
        else:
            raise ValueError(f"MHC type of {mhc} not recognised.")

        (self.tcr_com, self.mhc_com, self.tcr_VA_com, self.tcr_VB_com) = (
            self.mhc_tcr_com_calculator.calculate_centres_of_mass(
                tcr, save_aligned_as=save_aligned_as
            )
        )

        self.tcr_vector = self.get_tcr_vector(
            self.tcr_com, self.tcr_VA_com, self.tcr_VB_com
        )
        self.polarity = self.get_polarity(self.tcr_vector)
        self.scanning_angle, self.tcr_pitch_angle = self.get_tcr_docking_angles(
            self.tcr_vector, polarity_as_sign=polarity_as_sign
        )

    def __repr__(self):
        def polarity_to_str(polarity):
            return "canonical" if polarity == 0 else "reverse"

        return (
            f"TCR CoM: {self.tcr_com}\nMHC CoM: {self.mhc_com}\n"
            + f"TCR VA CoM: {self.tcr_VA_com}\nTCR VB CoM: {self.tcr_VB_com}\n"
            + f"Scanning angle: {np.degrees(self.scanning_angle)}\n"
            + f"Pitch angle: {np.degrees(self.tcr_pitch_angle)}\n"
            + f"Polarity: {polarity_to_str(self.polarity)}"
        )

    def to_dict(self):
        return {
            "tcr_com": self.tcr_com.tolist(),
            "mhc_com": self.mhc_com.tolist(),
            "tcr_VA_com": self.tcr_VA_com.tolist(),
            "tcr_VB_com": self.tcr_VB_com.tolist(),
            "scanning_angle": np.degrees(self.scanning_angle),
            "pitch_angle": np.degrees(self.tcr_pitch_angle),
            "polarity": self.polarity,
        }

    def get_tcr_vector(
        self, tcr_com: np.array, tcr_VA_com: np.array, tcr_VB_com: np.array
    ) -> np.array:
        """
        Calculates the vector from the VA centre of mass to the VB centre of mass,
        translated such that the vector originates at the total TCR's centre of mass and truncated to unit length.

        Args:
            tcr_com (np.array): Total TCR centre of mass
            tcr_VA_com (np.array): TCR VA centre of mass
            tcr_VB_com (np.array): TCR VB center of mass
        Returns:
            np.array: Unit vector from TCR CoM in direction of VA CoM to VB CoM
        """
        direction_vec = tcr_VB_com - tcr_VA_com
        tcr_vector = direction_vec / np.linalg.norm(direction_vec)
        return tcr_vector

    def get_tcr_docking_angles(
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
        scanning_angle = np.arccos(np.dot(xy_projection, np.asarray([0.0, 1.0])))
        phi = np.arccos(np.sqrt(1 - (tcr_vector[-1] ** 2)))
        if polarity_as_sign:
            scanning_angle = scanning_angle * ((-1) ** self.polarity)
        return scanning_angle, phi

    def get_polarity(self, tcr_vector: np.array) -> int:
        """
        Return the polarity of the TCR based on the TCR vector pointing from the VA to the VB CoM.
        If the x component is negative, ie the tcr_vector points from the alpha 2 chain to the
        alpha 1 chain of the MHC, the polarity is canonical (0). Otherwise the polarity is reverse (1).

        Args:
            tcr_vector (np.array): Unit vector pointing from VA CoM to VB CoM.

        Returns:
            int: 0 for canonical polarity, 1 for reverse polarity.
        """
        if tcr_vector[0] <= 0:
            return 0
        else:
            return 1
