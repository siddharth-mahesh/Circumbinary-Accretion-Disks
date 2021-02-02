import numpy as np
import matplotlib.pyplot as plt

gap_sizes = []

for i in range(5):
    res_file_name = "CoplanarLyapExpComputation"+str(i+1)+".txt"
    data = np.loadtxt(res_file_name)
    data_label = 'q = 0.'+str(i+1)
    r = data[:,0]
    max_lyap = np.array([max(data[j,1:]) for j in range(len(data))])
    time_scales = np.log10(1/max_lyap/2/np.pi)
    gap_sizes.append([(i+1)*0.1,r[np.abs(time_scales).argmin()]])
    plt.plot(r,time_scales,label = data_label)


plt.axhline(0,color = 'black')
plt.xlabel("Radial Distance")
plt.ylabel("Log(Instability Scale)")
plt.legend()
plt.savefig("InstabilityScale.png")
plt.show()

gap_sizes = np.array(gap_sizes)
np.savetxt("coplanar_gap_sizes.txt",gap_sizes)
plt.plot(gap_sizes[:,0],gap_sizes[:,1])
plt.ylabel("Gap Size")
plt.xlabel("Mass Ratio")
plt.savefig("GapSizes.png")
plt.show()


