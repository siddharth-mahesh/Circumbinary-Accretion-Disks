import numpy as np
import matplotlib.pyplot as plt

gap_sizes = []

for j in range(7):
    inc = (j)*15
    gap_sizes_for_this_inc = []
    for i in range(5):
        res_file_name = "InclinedLyapExpComputation_q"+str(i+1)+"_i"+str(inc)+".txt"
        print(res_file_name)
        data = np.loadtxt(res_file_name)
        data_label = 'q = 0.'+str(i+1)
        r = data[:,0]
        max_lyap = np.array([max(data[k,1:]) for k in range(len(data))])
        time_scales = np.log10(1/max_lyap/2/np.pi)
        gap_sizes_for_this_inc.append([(i+1)*0.1,r[np.abs(time_scales).argmin()]])
        plt.plot(r,time_scales,label = data_label)
    plt.axhline(0,color = 'black')
    plt.xlabel("Radial Distance")
    plt.ylabel("Log(Instability Scale)")
    plt.legend()
    figname = "InclinedInstabilityScale_"+"i"+str(inc)+".png"
    plt.savefig(figname)
    plt.show()
    gap_sizes.append(gap_sizes_for_this_inc)


gap_sizes = np.array(gap_sizes)
#print(gap_sizes)
#np.savetxt("inclined_gap_sizes.txt",gap_sizes)

print(gap_sizes)

for j in range(7):
    print("j = ",j)
    inc = (j)*15
    gaps_for_j = gap_sizes[j]
    print(gaps_for_j)
    print(gaps_for_j[:,0])
    data_label = "inc = "+str(inc)
    plt.plot(gaps_for_j[:,0],gaps_for_j[:,1],label = data_label)

plt.ylabel("Gap Size")
plt.xlabel("Mass Ratio")
plt.legend()
#plt.tight_layout()
plt.savefig("InclinedGapSizesByMassRatio.png")

plt.show()

inclinations = np.arange(0,105,15,dtype = int)

for j in range(5):
    print(gap_sizes[:,j,1])
    data_label = r"$\mu$ = 0."+str(j+1)
    plt.plot(inclinations,gap_sizes[:,j,1],label = data_label)

plt.ylabel("Gap Size")
plt.xlabel("Inclinations")
plt.legend()
#plt.tight_layout()
plt.savefig("InclinedGapSizesByInclination.png")



