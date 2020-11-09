import numpy as np
import matplotlib.pyplot as plt

for i in range(5):
    res_file_name = "LyapExpComputation"+str(i+1)+".txt"
    data = np.loadtxt(res_file_name)
    data_label = 'q = 0.'+str(i+1)
    r = data[:,0]
    max_lyap = np.array([max(data[j,1:]) for j in range(len(data))])
    plt.plot(r,-np.log10(max_lyap),label = data_label)


plt.axhline(np.log10(1./3e-4),label = "viscous scale")
plt.xlabel("Radial Distance")
plt.ylabel("Log(Instability Scale)")
plt.legend()
plt.savefig("InstabilityScale.png")
plt.show()
