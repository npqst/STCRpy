from ..tcr_processing import TCRParser

import numpy as np
import math
import sys


class TCRAngle(object):
    """
    Base class for calculating the orientation angles between the variable domains of a TCR.
    """

    dat_path = os.path.join(os.path.split(__file__)[0], "dat")
    sparser = TCRParser(QUIET=True)

    def _normalise(self, l):
        """
        Normalise a vector
        """
        a = np.array(l)
        return a / np.linalg.norm(a)

    def get_angles(self, tcr):
        return self.calculate_angles(tcr)


class abTCRAngle(TCRAngle):
    """
    Class to calculate the orientation angles between the variable domains of an abTCR.
    """

    def __init__(self):
        self._read_consensus()
        self._read_pc()
        self._read_coreset()

    def _read_consensus(self):
        """
        Read in the consensus structure coordinates.
        """
        self.Bconsensus = self.sparser.get_structure(
            "bconsensus", os.path.join(self.dat_path, "consensus_B.pdb")
        )
        self.Aconsensus = self.sparser.get_structure(
            "aconsensus", os.path.join(self.dat_path, "consensus_A.pdb")
        )

        self.Bcon_coords = {}
        for atom in self.Bconsensus.get_atoms():
            self.Bcon_coords[
                (atom.parent.parent.id,) + atom.parent.id[1:] + (atom.id,)
            ] = atom.get_coord()

        self.Acon_coords = {}
        for atom in self.Aconsensus.get_atoms():
            self.Acon_coords[
                (atom.parent.parent.id,) + atom.parent.id[1:] + (atom.id,)
            ] = atom.get_coord()

    def _read_pc(self):
        """
        Read in the plane vectors precalculated on the consensus structure
        """
        self.Apos = [
            list(map(float, x))
            for x in list(
                map(str.split, open(os.path.join(self.dat_path, "pcA.txt")).readlines())
            )
        ]
        self.Bpos = [
            list(map(float, x))
            for x in list(
                map(str.split, open(os.path.join(self.dat_path, "pcB.txt")).readlines())
            )
        ]

        # Use the antibody vector to allow for direct comparison
        self.cB = [
            -10 * 0.5 * self.Bpos[0][i] + 1 * 0.5 * self.Bpos[1][i] + self.Bpos[2][i]
            for i in range(3)
        ]
        self.cA = [
            6 * 0.5 * self.Apos[0][i] - 2 * 0.5 * self.Apos[1][i] + self.Apos[2][i]
            for i in range(3)
        ]

        # Define the plane vectors relative to "C"
        # On VA domain
        self.A1 = [self.cA[i] + self.Apos[0][i] for i in range(3)]
        self.A2 = [self.cA[i] + self.Apos[1][i] for i in range(3)]

        # On VB domain
        self.B1 = [self.cB[i] + self.Bpos[0][i] for i in range(3)]
        self.B2 = [self.cB[i] + self.Bpos[1][i] for i in range(3)]

    def _read_coreset(self):
        """
        Read in the coreset positions.
        """
        self.coresetA = [
            l.strip()[1:]
            for l in open(os.path.join(self.dat_path, "Acoreset.txt")).readlines()
        ]
        self.coresetB = [
            l.strip()[1:]
            for l in open(os.path.join(self.dat_path, "Bcoreset.txt")).readlines()
        ]

    def calculate_points(self, tcr):
        """
        Method to calculate the vector points for a TCR variable region.
        Note that this uses bio.pdb superimposer to do the structural fitting.
        @param tcr: A tcr structure object from TCRDB.TcrPDB
        """

        # grab the ca coordinates of the coreset positions of the beta and alpha domains.
        tcr_bcore, tcr_bcore_names = get_coords(
            tcr, [("residue", "B" + r) for r in self.coresetB], ["CA"]
        )
        tcr_acore, tcr_acore_names = get_coords(
            tcr, [("residue", "A" + r) for r in self.coresetA], ["CA"]
        )

        # make equivalent arrays for the consensus structures and the the variable domains.
        bcon_np, btcr_np = get_equivalent_arrays(
            self.Bcon_coords, list(self.Bcon_coords.keys()), tcr_bcore, tcr_bcore_names
        )
        acon_np, atcr_np = get_equivalent_arrays(
            self.Acon_coords, list(self.Acon_coords.keys()), tcr_acore, tcr_acore_names
        )

        B_rot, B_tran = superimpose(btcr_np, bcon_np).get_rotran()
        A_rot, A_tran = superimpose(atcr_np, acon_np).get_rotran()

        # Transform the coordinate system onto the variable region using the same rotation and translation
        Bpoints = [
            np.dot(self.cB, B_rot) + B_tran,
            np.dot(self.B1, B_rot) + B_tran,
            np.dot(self.B2, B_rot) + B_tran,
        ]
        Apoints = [
            np.dot(self.cA, A_rot) + A_tran,
            np.dot(self.A1, A_rot) + A_tran,
            np.dot(self.A2, A_rot) + A_tran,
        ]

        return Bpoints, Apoints

    def calculate_angles(self, tcr):
        """
        Method to calculate the orientation angles for a TCR.
        Note that this uses bio.pdb superimposer to do the structural fitting.
        @param tcr: A tcr structure object from TCRDB.TcrPDB
        """

        # grab the ca coordinates of the coreset positions of the heavy and light domains.
        tcr_bcore, tcr_bcore_names = get_coords(
            tcr, [("residue", "B" + r) for r in self.coresetB], ["CA"]
        )
        tcr_acore, tcr_acore_names = get_coords(
            tcr, [("residue", "A" + r) for r in self.coresetA], ["CA"]
        )

        # make equivalent arrays for the consensus structures and the the variable domains.
        bcon_np, btcr_np = get_equivalent_arrays(
            self.Bcon_coords, list(self.Bcon_coords.keys()), tcr_bcore, tcr_bcore_names
        )
        acon_np, atcr_np = get_equivalent_arrays(
            self.Acon_coords, list(self.Acon_coords.keys()), tcr_acore, tcr_acore_names
        )

        # superimpose the consensus domains onto the real structures.
        B_rot, B_tran = superimpose(btcr_np, bcon_np).get_rotran()
        A_rot, A_tran = superimpose(atcr_np, acon_np).get_rotran()

        # Transform the coordinate system onto the variable region using the same rotation and translation
        Bpoints = [
            np.dot(self.cB, B_rot) + B_tran,
            np.dot(self.B1, B_rot) + B_tran,
            np.dot(self.B2, B_rot) + B_tran,
        ]
        Apoints = [
            np.dot(self.cA, A_rot) + A_tran,
            np.dot(self.A1, A_rot) + A_tran,
            np.dot(self.A2, A_rot) + A_tran,
        ]

        # Create vectors with which to calculate angles between.
        C = self._normalise([Bpoints[0][i] - Apoints[0][i] for i in range(3)])
        Cminus = [-1 * x for x in C]
        A1 = self._normalise([Apoints[1][i] - Apoints[0][i] for i in range(3)])
        A2 = self._normalise([Apoints[2][i] - Apoints[0][i] for i in range(3)])
        B1 = self._normalise([Bpoints[1][i] - Bpoints[0][i] for i in range(3)])
        B2 = self._normalise([Bpoints[2][i] - Bpoints[0][i] for i in range(3)])
        dc = (
            sum([x**2 for x in [Bpoints[0][i] - Apoints[0][i] for i in range(3)]])
            ** 0.5
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
                    ["BA", "BC1", "AC1", "BC2", "AC2", "dc"],
                    [BA, BC1, AC1, BC2, AC2, dc],
                )
            )
        )


class gdangle(trangle):
    """
    Class to calculate the orientation angles between the variable domains of a gdTCR.
    """

    def __init__(self):
        self._read_consensus()
        self._read_pc()
        self._read_coreset()

    def _read_consensus(self):
        """
        Read in the consensus structure coordinates.
        """
        self.Dconsensus = self.sparser.get_structure(
            "dconsensus", os.path.join(self.dat_path, "consensus_D.pdb")
        )
        self.Gconsensus = self.sparser.get_structure(
            "gconsensus", os.path.join(self.dat_path, "consensus_G.pdb")
        )

        self.Dcon_coords = {}
        for atom in self.Dconsensus.get_atoms():
            self.Dcon_coords[
                (atom.parent.parent.id,) + atom.parent.id[1:] + (atom.id,)
            ] = atom.get_coord()
        self.Gcon_coords = {}
        for atom in self.Gconsensus.get_atoms():
            self.Gcon_coords[
                (atom.parent.parent.id,) + atom.parent.id[1:] + (atom.id,)
            ] = atom.get_coord()

    def _read_pc(self):
        """
        Read in the plane vectors precalculated on the consensus structure
        """
        self.Dpos = [
            list(map(float, x))
            for x in list(
                map(str.split, open(os.path.join(self.dat_path, "pcB.txt")).readlines())
            )
        ]
        self.Gpos = [
            list(map(float, x))
            for x in list(
                map(str.split, open(os.path.join(self.dat_path, "pcA.txt")).readlines())
            )
        ]

        # Use the antibody vector to allow for direct comparison
        self.cD = [
            -10 * 0.5 * self.Gpos[0][i] + 1 * 0.5 * self.Gpos[1][i] + self.Gpos[2][i]
            for i in range(3)
        ]
        self.cG = [
            6 * 0.5 * self.Dpos[0][i] - 2 * 0.5 * self.Dpos[1][i] + self.Dpos[2][i]
            for i in range(3)
        ]

        # Define the plane vectors relative to "C"
        # On VA domain
        self.G1 = [self.cG[i] + self.Gpos[0][i] for i in range(3)]
        self.G2 = [self.cG[i] + self.Gpos[1][i] for i in range(3)]

        # On VB domain
        self.D1 = [self.cD[i] + self.Dpos[0][i] for i in range(3)]
        self.D2 = [self.cD[i] + self.Dpos[1][i] for i in range(3)]

    def _read_coreset(self):
        """
        Read in the coreset positions.
        """
        self.coresetG = [
            l.strip()[1:]
            for l in open(os.path.join(self.dat_path, "Acoreset.txt")).readlines()
        ]
        self.coresetD = [
            l.strip()[1:]
            for l in open(os.path.join(self.dat_path, "Bcoreset.txt")).readlines()
        ]

    def calculate_points(self, tcr):
        """
        Method to calculate the vector points for a TCR variable region.
        Note that this uses bio.pdb superimposer to do the structural fitting.
        @param tcr: A tcr structure object from TCRDB.TcrPDB
        """

        # grab the ca coordinates of the coreset positions of the beta and alpha domains.
        tcr_gcore, tcr_gcore_names = get_coords(
            tcr, [("residue", "G" + r) for r in self.coresetB], ["CA"]
        )
        tcr_dcore, tcr_dcore_names = get_coords(
            tcr, [("residue", "D" + r) for r in self.coresetA], ["CA"]
        )

        # make equivalent arrays for the consensus structures and the the variable domains.
        gcon_np, btcr_np = get_equivalent_arrays(
            self.Gcon_coords, list(self.Gcon_coords.keys()), tcr_gcore, tcr_gcore_names
        )
        dcon_np, dtcr_np = get_equivalent_arrays(
            self.Dcon_coords, list(self.Dcon_coords.keys()), tcr_dcore, tcr_dcore_names
        )

        G_rot, G_tran = superimpose(btcr_np, gcon_np).get_rotran()
        D_rot, D_tran = superimpose(dtcr_np, dcon_np).get_rotran()

        # Transform the coordinate system onto the variable region using the same rotation and translation
        Gpoints = [
            np.dot(self.cG, G_rot) + G_tran,
            np.dot(self.G1, G_rot) + G_tran,
            np.dot(self.G2, G_rot) + G_tran,
        ]
        Dpoints = [
            np.dot(self.cD, D_rot) + D_tran,
            np.dot(self.D1, D_rot) + D_tran,
            np.dot(self.D2, D_rot) + D_tran,
        ]

        return Dpoints, Gpoints

    def calculate_angles(self, tcr):
        """
        Method to calculate the orientation angles for a TCR.
        Note that this uses bio.pdb superimposer to do the structural fitting.
        @param tcr: tcr structure object from TCRDB.TcrPDB
        """

        # grab the ca coordinates of the coreset positions of the heavy and light domains.
        tcr_gcore, tcr_gcore_names = get_coords(
            tcr, [("residue", "G" + r) for r in self.coresetG], ["CA"]
        )
        tcr_dcore, tcr_dcore_names = get_coords(
            tcr, [("residue", "D" + r) for r in self.coresetD], ["CA"]
        )

        # make equivalent arrays for the consensus structures and the the variable domains.
        gcon_np, btcr_np = get_equivalent_arrays(
            self.Gcon_coords, list(self.Gcon_coords.keys()), tcr_gcore, tcr_gcore_names
        )
        dcon_np, dtcr_np = get_equivalent_arrays(
            self.Dcon_coords, list(self.Dcon_coords.keys()), tcr_dcore, tcr_dcore_names
        )

        # superimpose the consensus domains onto the real structures.
        G_rot, G_tran = superimpose(btcr_np, gcon_np).get_rotran()
        D_rot, D_tran = superimpose(dtcr_np, dcon_np).get_rotran()

        # Transform the coordinate system onto the variable region using the same rotation and translation
        Gpoints = [
            np.dot(self.cG, G_rot) + G_tran,
            np.dot(self.G1, G_rot) + G_tran,
            np.dot(self.G2, G_rot) + G_tran,
        ]
        Dpoints = [
            np.dot(self.cD, D_rot) + D_tran,
            np.dot(self.D1, D_rot) + D_tran,
            np.dot(self.D2, D_rot) + D_tran,
        ]

        # Create vectors with which to calculate angles between.
        C = self._normalise([Gpoints[0][i] - Dpoints[0][i] for i in range(3)])
        Cminus = [-1 * x for x in C]
        D1 = self._normalise([Dpoints[1][i] - Dpoints[0][i] for i in range(3)])
        D2 = self._normalise([Dpoints[2][i] - Dpoints[0][i] for i in range(3)])
        G1 = self._normalise([Gpoints[1][i] - Gpoints[0][i] for i in range(3)])
        G2 = self._normalise([Gpoints[2][i] - Gpoints[0][i] for i in range(3)])
        dc = (
            sum([x**2 for x in [Gpoints[0][i] - Dpoints[0][i] for i in range(3)]])
            ** 0.5
        )

        # Projection of the D1 and G1 vectors onto the plane perpendicular to C.
        n_x = np.cross(D1, C)
        n_y = np.cross(C, n_x)

        tmpD_ = self._normalise([0, np.dot(D1, n_x), np.dot(D1, n_y)])
        tmpG_ = self._normalise([0, np.dot(G1, n_x), np.dot(G1, n_y)])

        # GD is the angle between the D1 and G1 vectors looking down C.
        DG = math.acos(np.dot(tmpG_, tmpD_))
        DG = DG * (180.0 / math.pi)

        # Find direction by computing cross products
        if np.dot(np.cross(tmpD_, tmpG_), [1, 0, 0]) < 0:
            DG = -DG

        # DC1 angle is the angle between the D1 and C vectors
        DC1 = math.acos(np.dot(D1, C))
        DC1 = DC1 * (180.0 / math.pi)

        # GC1 angle is the angle between the G1 and C vectors
        GC1 = math.acos(np.dot(G1, Cminus))
        GC1 = GC1 * (180.0 / math.pi)

        # DC2 angle is the angle between the D2 and C vectors
        DC2 = math.acos(np.dot(D2, C))
        DC2 = DC2 * (180.0 / math.pi)

        # GC2 angle is the angle between the G2 and C vectors
        GC2 = math.acos(np.dot(G2, Cminus))
        GC2 = GC2 * (180.0 / math.pi)

        # Return the angles and the separation distance.
        return dict(
            list(
                zip(
                    ["DG", "DC1", "GC1", "DC2", "GC2", "dc"],
                    [DG, DC1, GC1, DC2, GC2, dc],
                )
            )
        )
