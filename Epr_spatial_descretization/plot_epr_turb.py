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



##-------------------------------##
##------   varying dx -----------##
##-------------------------------##
dft=pd.read_csv('Turb_epr.csv')
dft.head()

dft1=pd.read_csv('Turb_epr_dx=0.25.csv')
dft1.head()


k1_t=0.1-dft['k1']
epr_t=dft['epr']
epr_tv=dft['epr_v']


k1_t_1=0.1-dft1['k1']
epr_t_1=dft1['epr']
epr_tv_1=dft1['epr_v']




fig, ax = plt.subplots()

plt.errorbar(k1_t_1,epr_t_1,yerr=np.sqrt(epr_tv_1), fmt='o',ecolor='b',ms='8',mec='b',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$\Delta=0.25$")
plt.plot(k1_t_1,epr_t_1,'--',c='b',linewidth=2)

plt.errorbar(k1_t,epr_t,yerr=np.sqrt(epr_tv), fmt='o',ecolor='g',ms='8',mec='g',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$\Delta=0.5$")
plt.plot(k1_t,epr_t,'--',c='g',linewidth=2)


plt.errorbar(k1_t_1,epr_t_1/4,yerr=np.sqrt(epr_tv_1/4), fmt='o',ecolor='r',ms='8',mec='r',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$\Delta=0.25, \frac{EPR}{4}$")
plt.plot(k1_t_1,epr_t_1/4,'--',c='r',linewidth=2)


plt.xlabel(r"$k-k^{'}$",fontsize=25)
plt.ylabel(r"$EPR$",fontsize=25)#plt.ylabel(r"$(k_{B}T) \times S$",fontsize=18)
plt.text(20.00,3e-2,r"$k=0.1$",fontsize=22)

plt.xlim(8e-2,115)

plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

ax.set_ylim(1e-5,5e2)
plt.xscale('log')
plt.yscale('log')
ax.set_yticks([1e-4,1e-2,1e0,1e2])
ax.set_yticklabels([r"$10^{-4}$",r"$10^{-2}$",r"$10^{0}$",r"$10^{2}$"],fontsize=25)
ax.tick_params(axis='both',which='major', direction='in',length=7.0)
ax.tick_params(axis='both',which='minor', direction='in',length=4.0)

plt.legend(loc='lower right',fontsize=20,frameon=True)


plt.savefig('epr_turb_dx.png', bbox_inches='tight',dpi=600)
print('done')



##-------------------------------##
##------   varying L  -----------##
##-------------------------------##
dft1=pd.read_csv('Turb_epr_L128.csv')
dft1.head()



k1_t_1=0.1-dft1['k1']
epr_t_1=dft1['epr']
epr_tv_1=dft1['epr_v']




fig, ax = plt.subplots()



plt.errorbar(k1_t,epr_t,yerr=np.sqrt(epr_tv), fmt='o',ecolor='g',ms='8',mec='g',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$L=64$")
plt.plot(k1_t,epr_t,'--',c='g',linewidth=2)

plt.errorbar(k1_t_1,epr_t_1,yerr=np.sqrt(epr_tv_1), fmt='o',ecolor='m',ms='8',mec='m',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$L=128$")
plt.plot(k1_t_1,epr_t_1,'--',c='m',linewidth=2)

plt.errorbar(k1_t_1,epr_t_1/4,yerr=np.sqrt(epr_tv_1/4), fmt='o',ecolor='brown',ms='8',mec='brown',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$L=128, \frac{EPR}{4}$")
plt.plot(k1_t_1,epr_t_1/4,'--',c='brown',linewidth=2)


plt.xlabel(r"$k-k^{'}$",fontsize=25)
plt.ylabel(r"$EPR$",fontsize=25)#plt.ylabel(r"$(k_{B}T) \times S$",fontsize=18)
plt.text(20.00,3e-2,r"$k=0.1$",fontsize=22)

plt.xlim(8e-2,115)

plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

ax.set_ylim(1e-5,5e2)
plt.xscale('log')
plt.yscale('log')
ax.set_yticks([1e-4,1e-2,1e0,1e2])
ax.set_yticklabels([r"$10^{-4}$",r"$10^{-2}$",r"$10^{0}$",r"$10^{2}$"],fontsize=25)
ax.tick_params(axis='both',which='major', direction='in',length=7.0)
ax.tick_params(axis='both',which='minor', direction='in',length=4.0)

plt.legend(loc='lower right',fontsize=20,frameon=True)


plt.savefig('epr_turbulent_L.png', bbox_inches='tight',dpi=600)
print('done')

plt.show()
