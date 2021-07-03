import numpy as np
import sys,os
import matplotlib.pyplot as plt

## edit this path to where you are storing the repo
repo_path = os.path.join(r"C:\Users","siddh","Documents\GitHub\Circumbinary-Accretion-Disks")

## no need to edit these paths
data_path = os.path.join(repo_path,"Data")

gap_sizes = []

for i in range(5):
    res_file_name = "InclinedLyapExpComputation_q"+str(i+1)+"_i0.txt"
    data = np.loadtxt(res_file_name)
    data_label = r'$\mu$ = 0.'+str(i+1)
    r = data[:,0]
    max_lyap = np.array([max(data[j,1:]) for j in range(len(data))])
    time_scales = np.log10(1/max_lyap/2/np.pi)
    #gap_sizes.append([(i+1)*0.1,r[np.abs(time_scales).argmin()]])
    plt.plot(r,time_scales,label = data_label)


plt.axhline(0,color = 'black',label = r'$\delta = 1$')
plt.xlabel("Radial Distance")
plt.ylabel(r'Log(Instability Scale ($\delta$))')
plt.legend()
plt.savefig("InstabilityScale.png")
plt.show()



coplanar_mass_ratios = np.arange(0.1,0.51,0.01)
gap_sizes = np.zeros(len(coplanar_mass_ratios))
for i in range(len(coplanar_mass_ratios)):
    res_file_name = "CoplanarLyapExpComputation"+str(int(100*coplanar_mass_ratios[i]))+".txt"
    data = np.loadtxt(res_file_name)
    r = data[:,0]
    max_lyap = np.array([max(data[j,1:]) for j in range(len(data))])
    time_scales = np.log10(1/max_lyap/2/np.pi)
    gap_sizes[i] = r[np.abs(time_scales).argmin()]
    
np.savetxt("coplanar_gap_sizes.txt",gap_sizes)
plt.plot(coplanar_mass_ratios,gap_sizes)
plt.ylabel(r'$r_\mathrm{L}$')
plt.xlabel(r'$\mu$')
plt.savefig("GapSizesRL.png")
plt.show()


numdatacoplanar = np.loadtxt(os.path.join(data_path,'aint0.dat'))
hires_mass_ratios = np.arange(0.05,0.51,0.01)
r_T = np.loadtxt("FluidGapSizeCoplanarHiRes.txt")

plt.plot(coplanar_mass_ratios,gap_sizes,color = 'black', label = r'$r_\mathrm{L}$')
plt.plot(hires_mass_ratios[5:],r_T[5:]*gap_sizes[-1]/r_T[-1], color = 'red',label = r'$\bar{r}_\mathrm{T}$')
plt.scatter(1/(1+1/numdatacoplanar[:,1]),gap_sizes[-1]*numdatacoplanar[:,4]/numdatacoplanar[0,4],edgecolors = 'black',facecolors = 'none',label = r'$\bar{r}_\mathrm{1\%}$')
plt.scatter(1/(1+1/numdatacoplanar[:,1]),gap_sizes[-1]*numdatacoplanar[:,3]/numdatacoplanar[0,3],edgecolors = 'red',facecolors = 'none',label = r'$\bar{r}_\mathrm{g}$')
plt.ylabel(r'$r_\mathrm{gap}$')
plt.xlabel(r'$\mu$')
plt.legend()
plt.savefig("GapSizesCoplanarRLR1p.png")
plt.show()


