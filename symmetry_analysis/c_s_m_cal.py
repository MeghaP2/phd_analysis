import csm
from csm.molecule.atom import Atom
from csm.calculations.data_classes import Operation
from csm.calculations import Exact
from openbabel import openbabel
import os
import MDAnalysis as mda
import numpy as np

# Load the universe from the trajectory file
#u = mda.Universe('test5.pdb', in_memory=True)
u = mda.Universe('solute.pdb','../first_1ms.xtc', in_memory=True)
#print("read the trajectory")
#u2 = u.select_atoms('resname BTA')
selections = [u.select_atoms(f"resid {i} and name C1 C2 C3 C4 C5 C6 C7 C8 C9 N1 N2 N3 O1 O2 O3") for i in range(1, 9)]
for ts in u.trajectory[100001:200000]:
    print(f"Frame {ts.frame}")  # Print current frame
    for j in range(len(selections)):
        #u2 = u.select_atoms(f"resid {i} and name C1 C2 C3 C4 C5 C6 C7 C8 C9 N1 N2 N3 Se1 Se2 Se3")
        pdb_filename = f"frame_{ts.frame}.pdb"
        
        try:
            selections[j].atoms.write(pdb_filename)
            reader = csm.molecule.molecule.MoleculeReader()
            mol = reader.from_file(pdb_filename)

            op = Operation("c3")
            calculation = Exact(op, mol)
            csm_value = calculation.calculate()

            # Write the result to file
            with open("from_100000_to_200001_frames_csm_values.txt", "a") as f:
                f.write(str(csm_value))
                f.write("\n")
        
        except Exception as e:
            print(f"Error in frame {ts.frame}: {e}")
        
        finally:
            # Ensure file deletion regardless of success or failure
            if os.path.exists(pdb_filename):
                os.remove(pdb_filename)
                print(f"Deleted temporary file: {pdb_filename}")

