
import warnings
import numpy as np

from ..tcr_processing import abTCR, MHCchain
from .TCRCoM import MHCI_TCRCoM, MHCII_TCRCoM


class TCRGeom():
    def __init__(self, tcr, save_aligned_as=False):
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
            warnings.warn('No MHC associated with TCR. Geometry cannot be calculated.')
        if not isinstance(tcr, abTCR):
            raise NotImplementedError(f'TCR MHC geometry only implemented for abTCR types, not {type(tcr)}')
                
        if isinstance(tcr.get_MHC()[0], MHCchain) or tcr.get_MHC()[0].MHC_type == 'MH1':
            self.mhc_tcr_com_calculator = MHCI_TCRCoM()
        elif tcr.get_MHC()[0].MHC_type == 'MH2':
            self.mhc_tcr_com_calculator = MHCII_TCRCoM()
        else:
            raise ValueError(f'MHC type of {tcr.get_MHC()} not recognised.')

        (self.tcr_com,
         self.mhc_com,
         self.tcr_VA_com,
         self.tcr_VB_com
         ) = self.mhc_tcr_com_calculator.calculate_centres_of_mass(tcr, save_aligned_as=save_aligned_as)

        self.tcr_vector = self.get_tcr_vector(self.tcr_com, self.tcr_VA_com, self.tcr_VB_com)
        self.scanning_angle, self.tcr_pitch_angle = self.get_tcr_docking_angles(self.tcr_vector, self.tcr_com)

    def __repr__(self):
        return (
            f'TCR CoM: {self.tcr_com}\nMHC CoM: {self.mhc_com}\nTCR VA CoM: {self.tcr_VA_com}\nTCR VB CoM: {self.tcr_VB_com}\nScanning angle: {np.degrees(self.scanning_angle)}\n Pitch angle: {np.degrees(self.tcr_pitch_angle)}'
            )
    
    def get_tcr_vector(self, tcr_com: np.array, tcr_VA_com: np.array, tcr_VB_com: np.array) -> np.array:
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

    def get_tcr_docking_angles(self, tcr_vector: np.array, tcr_com: np.array) -> tuple[np.array]:
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
            tcr_com (np.array): TCR centre of mass

        Returns:
            tuple[np.array]: Tuple containing the scanning angle and pitch angle of the TCR to the MHC
        """
        xy_projection = self.tcr_vector[:2] / np.linalg.norm(self.tcr_vector[:2])
        scanning_angle = np.arccos(np.dot(xy_projection, np.asarray([0., 1.])))
        phi = np.arccos(np.sqrt(1 - (tcr_vector[-1]**2)))
        return scanning_angle, phi
