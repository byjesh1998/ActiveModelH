import matplotlib.pyplot as plt
#from scipy.io import *
import pandas as pd
import numpy as np
from matplotlib import rc
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)




dfn=pd.read_csv('defect_slopes.csv')
dfn.head()


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


#laminar 
k1=0.1+dfn['k1']
lt=dfn['av_lt_p']


fig, ax = plt.subplots()

plt.xlabel(r"$k-k^{'}$",fontsize=25)
plt.ylabel(r"$Defect$ $lifetime$, $\tau_{L}$",fontsize=25)#plt.ylabel(r"$(k_{B}T) \times S$",fontsize=18)
plt.text(.9500,150,r"$k=0.1$",fontsize=20)
#plt.text(0.2500,0.4,r"$k^{'}=0.1$",fontsize=15)

#plt.text(.1,13,r"slope = %1.3f" %coef[0],fontsize=15)


#plt.xlim(0.050,1.2)
#plt.xlim(8e-2,15)

#ax.xaxis.set_major_locator(MultipleLocator(.2))
#ax.xaxis.set_minor_locator(MultipleLocator(.1))

#ax.yaxis.set_major_locator(MultipleLocator(2))
#ax.yaxis.set_minor_locator(MultipleLocator(1))

ax.tick_params(axis='both',which='both', direction='in')
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

plt.errorbar(k1,lt,yerr=0,fmt='s',ecolor='r',ms='8',mec='r',mfc='w',markeredgewidth=2,linewidth=1,capsize=3,label='Turbulent state')
plt.plot(k1,lt,'--',c='r')
plt.yscale('log')
#plt.xscale('log')
ax.tick_params(axis='both',which='major', direction='in',length=7.0)
ax.tick_params(axis='both',which='minor', direction='in',length=4.0)
#plt.xscale('log')
#plt.xlim(0.1,0.7)
#plt.xlim(0,1.2)
plt.text(.1,350.,r"$(c)$",fontsize=30)




##################################################################################################
##################################################################################################

left, bottom, width, height = [0.37, 0.35, 0.5, 0.5]
#ax2 = fig.add_axes([left, bottom, width, height])


dfn=pd.read_csv('msd_defect_k1_.01_p.csv')
dfn.head()

lt=dfn['t']
w=dfn['lt']

av=0
for i in range(len(lt)):
	add=lt[i]*w[i]
	av=av+add
	
av_lt=av/(sum(w))

print('Average life time:',av_lt)

s=[]
t=[]

for i in range(len(lt)):
	if w[i]>0:
		s.append(w[i])
		t.append(lt[i])
		
#plt.hist(s,bins=20,density=True)
###ax2.plot(t,s/sum(s),label=r"$k^{'}=-0.01$",linewidth=2,zorder=1e3)
#plt.xlim(0,100)

dfn=pd.read_csv('msd_defect_k1_.08_p.csv')
dfn.head()


lt=dfn['t']
w=dfn['lt']

av=0
for i in range(len(lt)):
	add=lt[i]*w[i]
	av=av+add
	
av_lt=av/(sum(w))

print('Average life time:',av_lt)

s=[]
t=[]

for i in range(len(lt)):
	if w[i]>0:
		s.append(w[i])
		t.append(lt[i])
		
#plt.hist(s,bins=20,density=True)
###ax2.plot(t,s/sum(s),label=r"$k^{'}=-0.08$",linewidth=2,zorder=1e2)


dfn=pd.read_csv('msd_defect_k1_.1_p.csv')
dfn.head()


lt=dfn['t']
w=dfn['lt']

av=0
for i in range(len(lt)):
	add=lt[i]*w[i]
	av=av+add
	
av_lt=av/(sum(w))

print('Average life time:',av_lt)

s=[]
t=[]

for i in range(len(lt)):
	if w[i]>0:
		s.append(w[i])
		t.append(lt[i])
		
#plt.hist(s,bins=20,density=True)
###ax2.plot(t,s/sum(s),label=r"$k^{'}=-0.1$",linewidth=2,zorder=5)

dfn=pd.read_csv('msd_defect_k1_.3_p_5e3.csv')
dfn.head()


lt=dfn['t']
w=dfn['lt']

av=0
for i in range(len(lt)):
	add=lt[i]*w[i]
	av=av+add
	
av_lt=av/(sum(w))

print('Average life time:',av_lt)

s=[]
t=[]

for i in range(len(lt)):
	if w[i]>0:
		s.append(w[i])
		t.append(lt[i])
		
#plt.hist(s,bins=20,density=True)
###ax2.plot(t,s/sum(s),label=r"$k^{'}=-0.3$",linewidth=2,zorder=1)



dfn=pd.read_csv('msd_defect_k1_1_p_5e3.csv')
dfn.head()


lt=dfn['t']
w=dfn['lt']

av=0
for i in range(len(lt)):
	add=lt[i]*w[i]
	av=av+add
	
av_lt=av/(sum(w))

print('Average life time:',av_lt)

s=[]
t=[]

for i in range(len(lt)):
	if w[i]>0:
		s.append(w[i])
		t.append(lt[i])
		
#plt.hist(s,bins=20,density=True)
###ax2.plot(t,s/sum(s),label=r"$k^{'}=-1$",linewidth=2,zorder=0)


#plt.xlim(.001,1050)
#plt.ylim(0,.5)
#plt.yscale('log')
#plt.xscale('log')
#plt.legend()

###ax2.set_ylabel(r"$P(\tau_{f})$")
###ax2.set_xlabel(r"$\tau_{f}$")











plt.savefig('defect_lt_k1.png', bbox_inches='tight',dpi=600)
plt.show()

