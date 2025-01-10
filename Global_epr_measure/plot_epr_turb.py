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
#       Varying k values         ##
##-------------------------------##
dft=pd.read_csv('Turb_epr_k=0.1.csv')
dft.head()


k1_t=0.1-dft['k1']
epr_t=dft['epr']
epr_tv=dft['epr_v']

dft2=pd.read_csv('Turb_epr_k=0.2.csv')
dft2.head()


k1_t2=0.2-dft2['k1']
epr_t2=dft2['epr']
epr_tv2=dft2['epr_v']



fig, ax = plt.subplots()
plt.errorbar(k1_t,epr_t,yerr=np.sqrt(epr_tv), fmt='o',ecolor='g',ms='8',mec='g',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$k=0.1$")
plt.plot(k1_t,epr_t,'--',c='g',linewidth=2)

plt.errorbar(k1_t2,epr_t2,yerr=np.sqrt(epr_tv2), fmt='o',ecolor='r',ms='8',mec='r',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$k=0.2$")
plt.plot(k1_t2,epr_t2,'--',c='r',linewidth=2)


#plot properties
plt.xlabel(r"$k-k^{'}$",fontsize=25)
plt.ylabel(r"$EPR$",fontsize=25)

plt.text(20.00,1e-2,r"$b=0.1$",fontsize=22)
plt.text(20.00,6e-2,r"$\eta=1.67$",fontsize=22)
plt.text(0.1,175.,r"$(a)$",fontsize=27)



plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

ax.set_ylim(1e-5,1e2)
plt.xscale('log')
plt.yscale('log')
ax.set_yticks([1e-4,1e-2,1e0,1e2])
ax.set_yticklabels([r"$10^{-4}$",r"$10^{-2}$",r"$10^{0}$",r"$10^{2}$"],fontsize=25)
ax.tick_params(axis='both',which='major', direction='in',length=7.0)
ax.tick_params(axis='both',which='minor', direction='in',length=4.0)

plt.legend(loc='lower right',fontsize=22,frameon=True)


plt.savefig('epr_turbulent.png', bbox_inches='tight',dpi=600)


##-------------------------------##
#       Varying b values         ##
##-------------------------------##

dfb0=pd.read_csv('Turb_epr_less.csv')
dfb0.head()


k1_b0=0.1-dfb0['k1']
epr_b0=dfb0['epr']
epr_bv0=dfb0['epr_v']

dfb=pd.read_csv('Turb_epr_b=0.2.csv')
dfb.head()


k1_b=0.1-dfb['k1'] 
epr_b=dfb['epr']
epr_bv=dfb['epr_v']

dfb2=pd.read_csv('Turb_epr_b=0.3.csv')
dfb2.head()


k1_b2=0.1-dfb2['k1']
epr_b2=dfb2['epr']
epr_bv2=dfb2['epr_v']



fig, ax = plt.subplots()
plt.errorbar(k1_b0[1:],epr_b0[1:],yerr=np.sqrt(epr_bv0[1:]), fmt='s',ecolor='g',ms='8',mec='g',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$b=0.1$")
plt.plot(k1_b0[1:],epr_b0[1:],'--',c='g',linewidth=2)

plt.errorbar(k1_b,epr_b ,yerr=np.sqrt(epr_bv), fmt='s',ecolor='m',ms='8',mec='m',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$b=0.2$")
plt.plot(k1_b,epr_b,'--',c='m',linewidth=2)

plt.errorbar(k1_b2,epr_b2,yerr=np.sqrt(epr_bv2), fmt='s',ecolor='b',ms='8',mec='b',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$b=0.3$")
plt.plot(k1_b2,epr_b2,'--',c='b',linewidth=2)


#plot properties
plt.xlabel(r"$k-k^{'}$",fontsize=25)
plt.ylabel(r"$EPR$",fontsize=25)

plt.text(20.00,2e-2,r"$k=0.1$",fontsize=22)
plt.text(20.00,9e-2,r"$\eta=1.67$",fontsize=22)
plt.text(0.1,175.,r"$(b)$",fontsize=27)



plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

ax.set_ylim(1e-5,1e2)
plt.xscale('log')
plt.yscale('log')
ax.set_yticks([1e-4,1e-2,1e0,1e2])
ax.set_yticklabels([r"$10^{-4}$",r"$10^{-2}$",r"$10^{0}$",r"$10^{2}$"],fontsize=25)
ax.tick_params(axis='both',which='major', direction='in',length=7.0)
ax.tick_params(axis='both',which='minor', direction='in',length=4.0)

plt.legend(loc='lower right',fontsize=22,frameon=True)


plt.savefig('epr_turbulent_bs.png', bbox_inches='tight',dpi=600)




##-------------------------------##
#       Varying eta values         ##
##-------------------------------##



dfe=pd.read_csv('Turb_epr_eta=3.csv')
dfe.head()


k1_e=0.1-dfe['k1'] 
epr_e=dfe['epr']
epr_ev=dfe['epr_v']

dfe2=pd.read_csv('Turb_epr_eta=5.csv')
dfe2.head()


k1_e2=0.1-dfe2['k1']
epr_e2=dfe2['epr']
epr_ev2=dfe2['epr_v']



fig, ax = plt.subplots()
plt.errorbar(k1_b0[1:],epr_b0[1:],yerr=np.sqrt(epr_bv0[1:]), fmt='^',ecolor='g',ms='8',mec='g',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$\eta=1.67$")
plt.plot(k1_b0[1:],epr_b0[1:],'--',c='g',linewidth=2)

plt.errorbar(k1_e,epr_e ,yerr=np.sqrt(epr_ev), fmt='^',ecolor='c',ms='8',mec='c',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$\eta=3$")
plt.plot(k1_e,epr_e,'--',c='c',linewidth=2)

plt.errorbar(k1_e2,epr_e2,yerr=np.sqrt(epr_ev2), fmt='^',ecolor='crimson',ms='8',mec='crimson',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label=r"$\eta=5$")
plt.plot(k1_e2,epr_e2,'--',c='crimson',linewidth=2)


#plot properties
plt.xlabel(r"$k-k^{'}$",fontsize=25)
plt.ylabel(r"$EPR$",fontsize=25)

plt.text(20.00,5e-2,r"$k=0.1$",fontsize=22)
plt.text(20.00,3e-1,r"$b=0.1$",fontsize=22)
plt.text(0.1,855.,r"$(c)$",fontsize=27)



plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

ax.set_ylim(1e-5,5e2)
plt.xscale('log')
plt.yscale('log')
ax.set_yticks([1e-4,1e-2,1e0,1e2])
ax.set_yticklabels([r"$10^{-4}$",r"$10^{-2}$",r"$10^{0}$",r"$10^{2}$"],fontsize=25)
ax.tick_params(axis='both',which='major', direction='in',length=7.0)
ax.tick_params(axis='both',which='minor', direction='in',length=4.0)

plt.legend(loc='lower right',fontsize=22,frameon=True)


plt.savefig('epr_turbulent_etas.png', bbox_inches='tight',dpi=600)

plt.show()
