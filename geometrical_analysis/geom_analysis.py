# step-1: import required modules
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PDB, GRO, XTC
import numpy
import numpy as np
from numpy import pi
from numpy.linalg import norm

# step-2: upload trajectory files

u1 = mda.Universe('solute.gro', 'solute.xtc')
u2 = u1.select_atoms('resname MOL') # solute resname

# step-3:

com_intra=[]
sd_intra=[]
ta_intra=[]
for ts in u1.trajectory:

    for i in np.arange(1, 3): # define residue ids
    

    
################ For Intramolecular com-com ####################
    
    #  Defining atom selections for NDI ring
        
    
        atom_sel_11 = u2.select_atoms('resid '+str(i)+' and name N2 C9 O2 C10 C11 C12 C18 C19 O3 N3 C28 O4 C29 C30 C31 C32 C33 C39 C40 O5').center_of_mass()
        atom_sel_12 = u2.select_atoms('resid '+str(i)+' and name N5 C45 O7 C46 C47 C53 C54 C55 O8 N6 C64 O9 C65 C66 C67 C68 C69 C70 C76 C76 O10').center_of_mass()

        diff_com = atom_sel_11 - atom_sel_12
        mag_dis_com = norm(diff_com)
        
        #Defining atom selections for Plane normal
    
        atoms_sel_1_N1 = u2.select_atoms('resid '+str(i)+' and name N2').positions
        atoms_sel_1_N2 = u2.select_atoms('resid '+str(i)+' and name N3').positions
        atoms_sel_1_N3 = u2.select_atoms('resid '+str(i)+' and name N5').positions
        atoms_sel_1_N4 = u2.select_atoms('resid '+str(i)+' and name N6').positions
    
        atoms_sel_1_C1 = u2.select_atoms('resid '+str(i)+' and name C11').positions
        atoms_sel_1_C2 = u2.select_atoms('resid '+str(i)+' and name C39').positions
        atoms_sel_1_C3 = u2.select_atoms('resid '+str(i)+' and name C53').positions
        atoms_sel_1_C4 = u2.select_atoms('resid '+str(i)+' and name C69').positions
        

        N21_1 = atoms_sel_1_N2 - atoms_sel_1_N1
        C21_1 = atoms_sel_1_C2 - atoms_sel_1_C1
        N21_2 = atoms_sel_1_N4 - atoms_sel_1_N3
        C21_2 = atoms_sel_1_C4 - atoms_sel_1_C3

    

        NDI_ring_vec_1 = np.cross(N21_1, C21_1)
        NDI_ring_mag_1 = norm(NDI_ring_vec_1)   
        
        NDI_ring_vec_2 = np.cross(N21_2, C21_2)
        NDI_ring_mag_2 = norm(NDI_ring_vec_2)
    
        n_cap_1 = NDI_ring_vec_1 / NDI_ring_mag_1
        flatten_n_cap_1 = n_cap_1.flatten()
        n_cap_2 = NDI_ring_vec_2 / NDI_ring_mag_2
        flatten_n_cap_2 = n_cap_2.flatten()
        
        dis_pi_pi = np.dot(diff_com, flatten_n_cap_1)
        mag_dis_pi_pi = norm(dis_pi_pi)
    
        com_intra.append(mag_dis_pi_pi) 
        

        
#################### For slip distance ################################
        per_to_dis_pi_pi = (dis_pi_pi)*(n_cap_1) #perpendicular to plane of NDI core of molecule
    
        slip_dis = diff_com - per_to_dis_pi_pi
    
        mag_slip_dis = norm(slip_dis)
    
        sd_intra.append(mag_slip_dis)
    

################ For twist angle #########################

        norm_N21_1 = norm(N21_1)
        norm_N21_2 = norm(N21_2)
    
        unit_vector_N21_1 = N21_1 / norm_N21_1
        flatten_unit_vector_N21_1 = unit_vector_N21_1.flatten()
        unit_vector_N21_2 = N21_2 / norm_N21_2
        flatten_unit_vector_N21_2 = unit_vector_N21_2.flatten()
    
        theta_value_vec = np.dot(flatten_unit_vector_N21_2, flatten_unit_vector_N21_1) 
        theta_value_mag = norm(theta_value_vec)
       
        twist_angle = (np.arccos(theta_value_mag))*(180/pi)
    
        ta_intra.append(twist_angle)

np.savetxt("sd_intra.txt", sd_intra, fmt="%1.2f")
np.savetxt("ta_intra.txt", ta_intra, fmt="%1.2f")
np.savetxt("com_intra.txt", com_intra, fmt="%1.2f")

