# Step-1 required modules

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PDB, GRO, XTC
import numpy as np
from numpy.linalg import norm

# Step-2 upload trajectories

u1 = mda.Universe('solute.gro','solute.xtc')
u2 = u1.select_atoms('resname MOL')


handedness = []
frame=[]
for ts in u1.trajectory:
    #print(ts.frame)
    for i in np.arange(1,2):
        # Define COMs for the NDI rings
        com1 = u2.select_atoms(
            f"resid {i} and name N2 C14 O2 C15 C16 C17 C20 C21 O4 N3 C30 O5 C31 C32 C33 C34 C35 C38 C39 O7 O3 O6"
        ).center_of_mass()
        com2 = u2.select_atoms(
            f"resid {i} and name N5 C53 O9 C54 C55 C56 C59 C60 O11 N6 C69 O12 C70 C71 C72 C73 C74 C77 C78 O14 O10 O13"
        ).center_of_mass()

        # Get positions of the atoms contributing to the transition vectors
        tvec1 = u2.select_atoms(f"resid {i} and name C15").positions[0]
        tvec2 = u2.select_atoms(f"resid {i} and name C31").positions[0]
        tvec3 = u2.select_atoms(f"resid {i} and name C54").positions[0]
        tvec4 = u2.select_atoms(f"resid {i} and name C70").positions[0]

        # Transition vectors
        vec1 = tvec2 - tvec1

        # Choose correct orientation of vec2 based on distances
        d1, d2 = norm(tvec1 - tvec3), norm(tvec1 - tvec4)
        d3, d4 = norm(tvec2 - tvec3), norm(tvec2 - tvec4)

        if d1 < d2 and d4 < d3:
            vec2 = tvec4 - tvec3
        elif d1 > d2 and d4 > d3:
            vec2 = tvec3 - tvec4
        else:
            continue

        # Cross product (μ1 × μ2)
        cp = np.cross(vec1, vec2)

        # Inter-chromophore vector
        r = com2 - com1

        # Normalized angle between dipoles
        norm_vec1 = vec1 / norm(vec1)
        norm_vec2 = vec2 / norm(vec2)
        cos_theta = np.dot(norm_vec1, norm_vec2)
        angle_mag = np.degrees(np.arccos(norm(cos_theta))) #, -1.0, 1.0)))

        # Scalar triple product
        mapping = np.dot(r, cp)

        # Signed handedness
        sign = np.sign(mapping)
        handedness.append(sign * angle_mag)
        frame.append(ts.frame)

np.savetxt("handedness_intra.txt", handedness, fmt="%1.2f")
