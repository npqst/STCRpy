from ..tcr_processing import TCRParser
from .reference_data import reference_data

import numpy as np
import math

from Bio.SVDSuperimposer import SVDSuperimposer


class TCRAngle:
    """
    Base class for calculating the orientation angles between the variable domains of a TCR.
    """

    def __init__(self, tcr=None):
        self.chain_types = self._get_TCR_chain_ids()
        self._get_reference_atoms()
        self._get_reference_vectors()

    def _get_TCR_chain_ids(self):
        raise NotImplementedError(
            "TCRAngle is a base class. Instantiate a descendent class, eg. abTCRAngle"
        )

    def _get_reference_atoms(self):
        parser = TCRParser.TCRParser(QUIET=True)
        reference_atoms = {}
        for chain_id in self.chain_types:
            reference_pdb_file = reference_data.reference_pdb_paths[chain_id]
            reference_structure = parser.get_tcr_structure(
                f"ref_{chain_id}", reference_pdb_file
            )
            reference_atoms[chain_id] = {
                atom.parent.id[1]: atom for atom in reference_structure.get_atoms()
            }
        self.reference_atoms = reference_atoms

    def _get_reference_vectors(self):
        """
        Read in the plane vectors precalculated on the reference structure
        """
        self.position_vectors = {
            chain_id: np.asarray(reference_data.reference_positions[chain_id])
            for chain_id in self.chain_types
        }

        self.plane_origin = {
            chain_id: np.dot(
                np.asarray(self.position_vectors[chain_id]).T,
                np.asarray(reference_data.reference_origin_scaling[chain_id]),
            )
            for chain_id in self.chain_types
        }

        # Define the plane vectors relative to "plane_origin"
        # On VA domain
        self.plane_1 = {
            chain_id: self.plane_origin[chain_id]
            + np.asarray(self.position_vectors[chain_id][0])
            for chain_id in self.chain_types
        }
        self.plane_2 = {
            chain_id: self.plane_origin[chain_id]
            + np.asarray(self.position_vectors[chain_id][1])
            for chain_id in self.chain_types
        }

    def _normalise(self, vec):
        """
        Normalise a vector
        """
        a = np.array(vec)
        return a / np.linalg.norm(a)

    def calculate_points(self, tcr, chain_id):
        """
        Method to calculate the vector points for a TCR variable region.
        Note that this uses bio.pdb superimposer to do the structural fitting.
        @param tcr: A tcr structure object from TCRDB.TcrPDB
        """

        def get_chain(tcr, chain_id):
            """Retrieve chain by chain type from TCR structure object

            Args:
                tcr (TCR structure object): TCR structure object
                chain_id (str): Chain identifier for TCR alpha, beta, gamma or delta chains.
                Should be in set {A, B, G, D}

            Returns:
                Chain structure object
            """
            if chain_id == "A":
                return tcr.get_VA()
            elif chain_id == "B":
                return tcr.get_VB()
            elif chain_id == "G":
                return tcr.get_VG()
            elif chain_id == "D":
                return tcr.get_VG()
            else:
                raise ValueError(f"Unrecognised TCR chain type {chain_id}")

        tcr_atoms = {
            res_id: get_chain(tcr, chain_id)[res_id]["CA"]
            for res_id in reference_data.reference_residues[chain_id]
            if res_id in get_chain(tcr, chain_id)
        }

        coordinates = np.asarray(
            [atom.get_coord() for res_id, atom in tcr_atoms.items()]
        )

        reference_coordinates = np.asarray(
            [
                atom.get_coord()
                for res_id, atom in self.reference_atoms[chain_id].items()
                if res_id in tcr_atoms
            ]
        )

        assert coordinates.shape == reference_coordinates.shape

        superimposer = SVDSuperimposer()

        def superimpose(coord, ref_coord):
            superimposer.set(coord, ref_coord)
            superimposer.run()
            return superimposer.get_rotran()

        rotation, translation = superimpose(coordinates, reference_coordinates)

        transformed_points = (
            np.dot(self.plane_origin[chain_id], rotation) + translation,
            np.dot(self.plane_1[chain_id], rotation) + translation,
            np.dot(self.plane_2[chain_id], rotation) + translation,
        )
        return transformed_points

    def calculate_angles(self, tcr):
        """
        Method to calculate the orientation angles for a TCR.
        Note that this uses bio.pdb superimposer to do the structural fitting.
        @param tcr: A tcr structure object from TCRDB.TcrPDB
        """

        plane_origin, plane_1, plane_2 = {}, {}, {}
        for chain_type in self.chain_types:
            plane_origin[chain_type], plane_1[chain_type], plane_2[chain_type] = (
                self.calculate_points(tcr, chain_type)
            )

        # Create vectors with which to calculate angles between.
        C = self._normalise(
            plane_origin[self.chain_types[1]] - plane_origin[self.chain_types[0]]
        )
        Cminus = [-1 * x for x in C]
        A1 = self._normalise(
            plane_1[self.chain_types[0]] - plane_origin[self.chain_types[0]]
        )
        A2 = self._normalise(
            plane_2[self.chain_types[0]] - plane_origin[self.chain_types[0]]
        )
        B1 = self._normalise(
            plane_1[self.chain_types[1]] - plane_origin[self.chain_types[1]]
        )
        B2 = self._normalise(
            plane_2[self.chain_types[1]] - plane_origin[self.chain_types[1]]
        )
        dc = np.sqrt(
            (
                (plane_origin[self.chain_types[1]] - plane_origin[self.chain_types[0]])
                ** 2
            ).sum()
        )

        # Projection of the A1 and B1 vectors onto the plane perpendicular to C.
        n_x = np.cross(A1, C)
        n_y = np.cross(C, n_x)

        tmpA_ = self._normalise([0, np.dot(A1, n_x), np.dot(A1, n_y)])
        tmpB_ = self._normalise([0, np.dot(B1, n_x), np.dot(B1, n_y)])

        # BA is the angle between the A1 and B1 vectors looking down C.
        BA = math.acos(np.dot(tmpA_, tmpB_))
        BA = BA * (180.0 / math.pi)

        # Find direction by computing cross products
        if np.dot(np.cross(tmpA_, tmpB_), [1, 0, 0]) < 0:
            BA = -BA

        # AC1 angle is the angle between the A1 and C vectors
        AC1 = math.acos(np.dot(A1, C))
        AC1 = AC1 * (180.0 / math.pi)

        # BC1 angle is the angle between the B1 and C vectors
        BC1 = math.acos(np.dot(B1, Cminus))
        BC1 = BC1 * (180.0 / math.pi)

        # AC2 angle is the angle between the A2 and C vectors
        AC2 = math.acos(np.dot(A2, C))
        AC2 = AC2 * (180.0 / math.pi)

        # BC2 angle is the angle between the B2 and C vectors
        BC2 = math.acos(np.dot(B2, Cminus))
        BC2 = BC2 * (180.0 / math.pi)

        # Return the angles and the separation distance.
        return dict(
            list(
                zip(
                    [
                        f"{self.chain_types[1]}{self.chain_types[0]}",
                        f"{self.chain_types[1]}C1",
                        f"{self.chain_types[0]}C1",
                        f"{self.chain_types[1]}C2",
                        f"{self.chain_types[0]}C2",
                        "dist",
                    ],
                    [BA, BC1, AC1, BC2, AC2, dc],
                )
            )
        )

    def get_angles(self, tcr):
        # handle cases where multiple TCRs are passed as one TCR structure model
        try:
            return {t.id: self.calculate_angles(t) for t in tcr.get_TCRs()}
        # else return angles for single TCR
        except AttributeError:
            return self.calculate_angles(tcr)


class abTCRAngle(TCRAngle):
    """
    Class to calculate the orientation angles between the variable domains of an abTCR.
    """

    def __init__(self):
        super().__init__()

    def _get_TCR_chain_ids(self):
        return ["A", "B"]


class gdTCRAngle(TCRAngle):
    """
    Class to calculate the orientation angles between the variable domains of a gdTCR.
    """

    def __init__(self):
        super().__init__()

    def _get_TCR_chain_ids(self):
        return ["G", "D"]


class dbTCRAngle(TCRAngle):
    """
    Class to calculate the orientation angles between the variable domains of a dbTCR.
    """

    def __init__(self):
        super().__init__()

    def _get_TCR_chain_ids(self):
        return ["D", "B"]
