import matplotlib.pyplot as plt
from scipy.io import *
import pandas as pd
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import rc




##define figure parameters
params = {'legend.fontsize': 6,
          'axes.linewidth': 1.5,
          'axes.labelsize': 6.5,
          'xtick.labelsize': 6,
          'ytick.labelsize': 6,	
          'text.usetex': True
          }
plt.rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})



dfl=pd.read_csv('epr_dx.csv')
dfl.head()


L=64

dx=dfl['dx']
epr=dfl['epr']
epr_v=dfl['epr_v']


fig, ax = plt.subplots()
plt.errorbar(dx,epr,yerr=np.sqrt(epr_v),fmt='s',ecolor='b',ms='7',mec='b',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label='Turbulent state')


plt.xlabel(r"$\Delta$",fontsize=25)
plt.ylabel(r"$EPR$",fontsize=25)#plt.ylabel(r"$(k_{B}T) \times S$",fontsize=18)
plt.text(0.2500,0.3,r"$k=0.1$",fontsize=22)
plt.text(0.2500,0.2,r"$k^{'}=0.1$",fontsize=22)

ax.set_ylim(0.1,0.9)
ax.set_yticks([0.2,0.4,0.6,0.8],fontsize=22)
ax.set_yticklabels([r"$0.2$",r"$0.4$",r"$0.6$",r"$0.8$"],fontsize=22)
#ax.set_xlim(0.2,0.6)
ax.set_xticks([0.25,0.3,0.35,0.4,0.45,0.5],fontsize=22)
ax.set_xticklabels([r"$0.25$",r"$0.3$",r"$0.35$",r"$0.4$",r"$0.45$",r"$0.5$"],fontsize=22)



left, bottom, width, height = [0.54, 0.5, 0.35, 0.35]
ax2 = fig.add_axes([left, bottom, width, height])
ax2.errorbar(dx,epr*dx**2/(L**2),yerr=np.sqrt(epr_v)*dx**2/(L**2),fmt='s',ecolor='b',ms='5',mec='b',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label='Turbulent state')
ax2.set_ylim(0.5e-5,2e-5)
ax2.set_yticks([1e-5,1.5e-5,2e-5])
ax2.set_yticklabels([r"$1$",r"$1.5$",r"$2$"],fontsize=20)


ax2.set_xticks([0.25,0.35,0.45])
ax2.set_xticklabels([r"$0.25$",r"$0.35$",r"$0.45$"],fontsize=20)

ax2.text(0.254,1.7e-5,r"$\times 10^{-5}$",fontsize=20)
ax2.set_ylabel(r"$EPR/N^{2}$",fontsize=22)
ax2.set_xlabel(r"$\Delta$",fontsize=22)

plt.savefig('epr_turb_dx2.png', bbox_inches='tight',dpi=600)
plt.show()
