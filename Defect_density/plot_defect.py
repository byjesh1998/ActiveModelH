import matplotlib.pyplot as plt
from scipy.io import *
import pandas as pd
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import rc
from numpy.polynomial import polynomial as P





dft=pd.read_csv('Turb_full_defects.csv')
dft.head()

# dft1=pd.read_csv('Turb_full_defects_4.csv')
# dft1.head()

#Turb_full_defects_dx=0.25


#turb
k1_t=0.1-dft['k1']
def_t=dft['def']
def_v_tv=dft['def_v']

# k1_t_1=0.1-dft1['k1']
# def_t_1=dft1['def']
# def_v_tv_1=dft1['def_v']


Lsq=64*64
#plt.scatter(k1_t,epr_t)

params = {'legend.fontsize': 6,
          #'text.latex.preamble': [r"\usepackage{amstext}",],
          #'font',**{'family':'sans-serif','sans-serif':['Helvetica']},
          'axes.linewidth': 1.5,
          'axes.labelsize': 6.5,
          #'text.fontsize': 4,
          'xtick.labelsize': 6,
          'ytick.labelsize': 6,	
          'text.usetex': True
          #'figure.figsize': fig_size
          }
plt.rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
fig, ax = plt.subplots()
#yr=np.sqrt(def_v_tv)/Lsq

plt.errorbar(k1_t,(def_t/4.096),yerr=np.sqrt(def_v_tv/(4.096)), fmt='s',ecolor='b',ms='6',mec='b',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label='Turbulent state')
#plt.errorbar(k1_t_1,(def_t_1/4.096),yerr=np.sqrt(def_v_tv_1/4.096), fmt='s',ecolor='b',ms='5',mec='b',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label='Turbulent state')

#plt.plot(k1_t_1,def_t_1/4.096,'r','--')

##fitting

x = np.log(k1_t[13:25])
y = np.log(def_t[13:25])

coef = np.polyfit(x,y,1)
print(coef)
poly1d_fn = np.poly1d(coef) 
# poly1d_fn is now a function which takes in x and returns an estimate for y
#ax2.set_title(r"$\langle{\Delta\mathbf{r}^2(t)\rangle} = %1.3f t^{%1.3f}$" %(exp(coef[1]),coef[0]))
#plt.plot( np.exp(x), np.exp(poly1d_fn(x)), '--k',linewidth=2)



#plot properties
plt.xlabel(r"$k-k^{'}$",fontsize=25)
plt.ylabel(r"$Defect$ $density$, $\rho_{D}$",fontsize=25)
plt.text(0.29,3,r"$k=0.1$",fontsize=25)

plt.text(0.1,125.,r"$(b)$",fontsize=27)
#plt.text(0.1,75.,r"$\times 10^{-3}$",fontsize=20)

ax.tick_params(axis='both',which='major', direction='in',length=7.0)
ax.tick_params(axis='both',which='minor', direction='in',length=4.0)


#plt.xlim(0.050,1.2)
#plt.ylim(0,60)

#ax.xaxis.set_major_locator(MultipleLocator(1))
#ax.xaxis.set_minor_locator(MultipleLocator(.5))

#ax.yaxis.set_major_locator(MultipleLocator(10))
#ax.yaxis.set_minor_locator(MultipleLocator(5))

#ax.tick_params(axis='both',which='both', direction='in')
plt.xticks(fontsize=25)
plt.yticks(fontsize=27)
#plt.legend(loc='upper left',fontsize=15)
#plt.plot(k1_t,k1_t**14)
plt.yscale('log')
plt.xscale('log')
ax.set_yticks([10,100])
ax.set_yticklabels([r"$10^{-2}$",r"$10^{-1}$"],fontsize=25)

left, bottom, width, height = [0.54, 0.2, 0.35, 0.35]
ax2 = fig.add_axes([left, bottom, width, height])
#ax2.plot(k1_t,(def_t),'r','--')
ax2.errorbar(k1_t,(def_t/4.096),yerr=np.sqrt(def_v_tv/4.096), fmt='s',ecolor='b',ms='6',mec='b',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label='Turbulent state')
#ax2.plot(k1_t,(def_t),'r','--')
ax2.set_xlim(0.05,1.1)
ax2.text(0.07,71.,r"$\times 10^{-3}$",fontsize=20)
#ax2.set_xticks(fontsize=18)
ax2.set_yticks([25,50,75])
ax2.set_yticklabels([r"$25$",r"$50$",r"$75$"],fontsize=20)

ax2.set_xticks([.25,.5,.75,1])
ax2.set_xticklabels([r"$0.25$",r"$0.5$",r"$0.75$",r"$1$"],fontsize=20)


plt.savefig('defect_turbulent.png', bbox_inches='tight',dpi=600)
print('done')
plt.show()
