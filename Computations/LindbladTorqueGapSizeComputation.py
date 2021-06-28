print("imports")
import sys,os
import numpy as np

print("setting up paths")
## edit only this path to where you are storing the repo
repo_path = os.path.join(r"C:\Users","siddh","Documents\GitHub\Circumbinary-Accretion-Disks")

## no need to edit these paths
formulae_path = os.path.join(repo_path,"Formulae")
sys.path.insert(0,formulae_path)
results_path = os.path.join(repo_path,"Results")

import ResonantTorquePicture as rtp

print("initializing arrays")
## set choices for mass ratios and radio

mass_ratios = np.arange(0.1,0.5+0.1,0.1)
inclinations = np.arange(0,1,0.01)*np.pi/2
eccentricities = np.arange(0.0,0.5,0.01)

## test case: R = 1e4

alpha = 1e-2
chi = 1e-1

## Inclined non eccentric cases

for i in range(len(mass_ratios)):
    print("\n q = %.2e \n"%mass_ratios[i])
    res_file_name = "FluidGapSizeByInclination_q"+str(i+1)+".txt"
    res_file = os.path.join(results_path,res_file_name)
    xgap = np.zeros(len(inclinations))
    for j in range(len(inclinations)):
        rgap = 0
        rgapparams = np.zeros(6)
        incl = inclinations[j]
        for m in range(5,0,-1):
            rgapguess = rtp.LR_location(m,1)
            params = [rgapguess,mass_ratios[i],incl,0.0,m,1]
            fluidparams = [params,alpha,chi]
            zeta = rtp.zeta_T(fluidparams)
            if zeta > 1:
                if rgapguess > rgap:
                    rgap = rgapguess
                    rgapparams = params
                    
        rgapguess = rtp.LR_location(2,2)
        params = [rgapguess,mass_ratios[i],incl,0.0,2,2]
        fluidparams = [params,alpha,chi]
        zeta = rtp.zeta_T(fluidparams)
        if zeta > 1:
            if rgapguess > rgap:
                rgap = rgapguess
                rgapparams = params
        if rgap != 0:
            xgap[j] = rtp.find_rgap( rgapparams , [rgap+1e-10,rgap+2])
        #xgap[j] = rgap
    np.savetxt(res_file,xgap)

## non-inclined eccentric case
    
    


