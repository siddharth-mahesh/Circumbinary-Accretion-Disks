import numpy as np
import matplotlib.pyplot as plt

mass_ratios = np.arange(0.1,0.6,0.1)
inclinations = np.arange(0,1,0.01)*90

for i in range(len(mass_ratios)):
    gapsizepath = "FluidGapSizeByInclination_q"+str(i+1)+".txt"
    plot_label = r'$\mu$ = 0.'+str(i+1)
    xgap = np.loadtxt(gapsizepath)
    plt.plot(inclinations,xgap,label = plot_label)

plt.xlabel(r'i($^\circ$)')
plt.ylabel(r'gap size ($a$)')
plt.legend()
plt.savefig('FluidGapSizeByInclination.png',dpi = 300)
plt.show()
    
