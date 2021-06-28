import numpy as np
import matplotlib.pyplot as plt

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
plt.ylabel("Gap Size")
plt.xlabel("Mass Ratio")
plt.savefig("GapSizes.png")
plt.show()


