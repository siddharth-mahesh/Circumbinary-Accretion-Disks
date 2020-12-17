print("imports")
import sys,os
import numpy as np

print("setting up paths")
## edit only this path to where you are storing the repo
repo_path = os.path.join(r"C:\Users","Siddharth Mahesh","Documents\GitHub\Circumbinary-Accretion-Disks")

## no need to edit these paths
formulae_path = os.path.join(repo_path,"Formulae")
sys.path.insert(0,formulae_path)
results_path = os.path.join(repo_path,"Results")

print("initializing arrays")
## set choices for mass ratios and radio

mass_ratios = np.arange(0.0,0.5+0.1,0.1)
avg_radii = np.arange(1.10,3,0.001)

print("primary computation")
## compute the jacobi constants

import StabilityMatrix as sm
from scipy.linalg import eigvals

mmax = 3
mmin = 1
res_index_file = open(os.path.join(results_path,"LyapExpComputationLabels.txt"),"w")
res_index_file.write("r_avg \t l0 \t l1 \t l2 \t l3")
res_index_file.close()

test_file_name = "test_keplerian.txt"
test_file = open(os.path.join(results_path,test_file_name),"w")

for i in range(len(mass_ratios)):
    q = mass_ratios[i]
    print("q = ",q)
    res_file_name = "LyapExpComputation"+str(int(10*q))+".txt"
    res_file = open(os.path.join(results_path,res_file_name),"w")
    print("mass ratio : %f"%(q))
    for j in range(len(avg_radii)):
        print("r = ", avg_radii[j])
        params = [avg_radii[j],mass_ratios[i],mmax,mmin]
        stability_mat = sm.K(params)
        #print(stability_mat)
        lyapexps = eigvals(stability_mat)
        #print(lyapexps)
        if i == 0:
            test_file.write("%f \t %f \t %f \t %f \t %f \n"%(avg_radii[j],lyapexps[0].imag,lyapexps[1].imag,lyapexps[2].imag,lyapexps[3].imag))
        res_file.write("%f \t %f \t %f \t %f \t %f \n"%(avg_radii[j],lyapexps[0].real,lyapexps[1].real,lyapexps[2].real,lyapexps[3].real))
    res_file.close()
test_file.close()
print("done")