import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt

from MDAnalysis.lib.distances import self_distance_array

u = mda.Universe("seq.gro", "whole.xtc")

ow = u.select_atoms("name OW1")

dist_overall = []
#l=int(u.dimensions[0])
#print(l)
for ts in u.trajectory:
    coords = ow.positions
    distances = self_distance_array(coords,box=u.dimensions)
    dist_overall.append(np.sort(distances)[:2000])

dist_overall = np.array(dist_overall)

print("distance calculation is done")

def inc_rdf_plot(sorted_dis):
    overall_val=[]
    for i in range(sorted_dis.shape[1]): # loop over 30 dis
        val=[]
        for j in range(sorted_dis.shape[0]):
            val.append(sorted_dis[j,i])
        overall_val.append(val)
        
    return overall_val

values=inc_rdf_plot(dist_overall)
values=np.array(values)

print(f"shape of values is: {values.shape}")

overall_x=[]
overall_hist=[]
for i in range(values.shape[0]):
    min_val=0 
    max_val=max(values[i,:])
    bin_size=0.01
    num_bins=np.arange(min_val, max_val+bin_size, bin_size)
    hist, edges = np.histogram(values[i], bins=num_bins)
    #r = 0.5 * (edges[1:] + edges[:-1])
    overall_hist.append(hist)
    overall_x.append(num_bins)

print("--------------------------------------")
print(f"(the shape of x is {np.array(overall_x).shape} and hist is {np.array(overall_hist).shape})")
print("--------------------------------------")

print("normalization part has started")

l=u.dimensions[0]
n_ow=len(ow)
rho=(n_ow)/l**3
overall_g_r=[]
num_frames=len(u.trajectory)

#--------------- plotting all RDFs separately --------------------------------
#for j in range(len(overall_hist)):
#    g_r=np.zeros(len(overall_x[j]))
#    #print(g_r)
#    for idx, k in enumerate(overall_hist[j]):
#        r=overall_x[j][idx]+bin_size/2
#        #print(bin_size)
#        norm=4*np.pi*(r**2)*(bin_size)*rho*num_frames
#        #print(norm)
#        g_r[idx]=k/norm
#    overall_g_r.append(g_r)
#-----------------------------------------------------------------------------
overall_x2=[]
k=0
j=0
while j < len(overall_hist):
    x_val=np.unique([x for sublist in overall_x[k:k+50] for x in sublist])
    #print(f"x_val are: {x_val}")
    g_r=np.zeros(len(x_val))
    for l in range(k,k+50):
        #k+=l      
        for idx, m in enumerate(overall_hist[l]):
            #print(f"overall_x[l][idx] is: {overall_x[l][idx]}")
            r=overall_x[l][idx]+bin_size/2
            #print(f"r is: {r}")
            #print(f"idx is: {idx}")
            norm=4*np.pi*(r[l][idx]**2)*(bin_size)*(rho)*num_frames/2
            g_r[idx]+=m/norm
    overall_g_r.append(g_r/50)
    overall_x2.append(x_val)
    j+=50
    k+=50
            

print("normalization is done")

plt.figure(figsize=(10,8))
colors = plt.cm.RdYlBu(np.linspace(0,1,len(overall_g_r)))
for i in range(len(overall_g_r)):
    plt.plot(overall_x2[i], overall_g_r[i],color=colors[i],alpha=1,linewidth=2, label=f"{i+1}")


# axis labels
plt.xlabel("r (Å)", fontsize=16)
plt.ylabel("O-O iRDF", fontsize=16)
plt.xlim(2.2,6)
# tick font sizes
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# grid
plt.grid(True, linestyle="--", alpha=0.8)

# legend outside plot
plt.legend(fontsize=8, bbox_to_anchor=(1.02,1), loc="upper left", ncol=4, borderaxespad=0)

# prevent cutoff
plt.tight_layout()
#plt.savefig("irdf.jpg", dpi=600)
plt.savefig("irdf2.jpg", dpi=600)

