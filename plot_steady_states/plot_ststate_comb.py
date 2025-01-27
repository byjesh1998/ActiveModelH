import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.ticker import  *
from matplotlib.pyplot import  *
from matplotlib import animation
from matplotlib import cm
from mpl_toolkits import mplot3d
from matplotlib import rc
import matplotlib as mpl



def vel_field(nx,ny,psi):
    dx=0.5
    epsilon = 1.0e-8

    kx = np.empty(nx)
    ky = np.empty(ny)
    
    kx[0:int(nx/2)] = np.float64(np.arange(0,int(nx/2))) *2*np.pi/(np.float64(nx)*dx)
    kx[int(nx/2):nx] = np.float64(np.arange(-int(nx/2),0))* 2*np.pi/(np.float64(nx)*dx)

    ky[0:ny] = kx[0:ny]
    
    kx[0] = epsilon
    ky[0] = epsilon
    kx, ky = np.meshgrid(kx, ky)
    
    psi_ft=np.fft.fft2(psi)
    

    vx_ft= 1j*ky[:,:]* psi_ft
    vy_ft= -1j*kx[:,:]* psi_ft


    vx=np.real(np.fft.ifft2(vx_ft))
    vy=np.real(np.fft.ifft2(vy_ft))
    	
    
    return vx,vy


fig_width_pt  = 116.
inches_per_pt = 1./72
fig_width 	  = fig_width_pt*inches_per_pt
fig_size	 	  = [3*fig_width, 2*fig_width]

params = {'legend.fontsize': 6,
          'axes.linewidth': 1.5,
          'axes.labelsize': 6.5,
          'xtick.labelsize': 6,
          'ytick.labelsize': 6,	
          'text.usetex': True,
          'figure.figsize': fig_size}
rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})



#plot 1
left   = 5.e-2
bottom = 5.50e-1 
width  = 2.8e-1
height = 4.0e-1

#plot 2
leftb   = 3.6e-1
bottomb = bottom
widthb  = width
heightb = height

#plot 3
leftc   = 6.68e-1
bottomc = bottom
widthc  = width
heightc = height

#plot 4
leftd   = 5e-2
bottomd = 0.83e-1
widthd  = 2.8e-1
heightd = 4.0e-1

#plot 5
lefte   = 3.6e-1
bottome = bottomd
widthe  = width
heighte = height

#plot 6
leftf   = 6.685e-1
bottomf = bottomd
widthf  = width
heightf = height


x = np.linspace(0.0, 128.0, 128)
y = np.linspace(0.0, 128.0, 128)
X, Y = np.meshgrid(x,y)


#laminar data
k1=-0.003
eta=1.67

phi_lam=np.loadtxt('phi_laminar.txt')

psi_lam=np.loadtxt('psi_laminar.txt')/abs(k1/eta)


k1=-0.0035
phi_vor=np.loadtxt('phi_vor.txt')

psi_vor=np.loadtxt('psi_vor.txt')/abs(k1/eta)

k1=-0.06
phi_turb=np.loadtxt('phi_turb.txt')

psi_turb=np.loadtxt('psi_turb.txt')/abs(k1/eta)



#plot 1
dum     = abs(phi_lam.flatten())
phi_max = max(dum)


fig = figure()
ax  = axes([left, bottom, width, height])
ax.set_xticks([])
ax.set_yticks([])
density = ax.pcolor(phi_lam, cmap='RdBu',vmin=-phi_max,vmax=phi_max)

props = dict(boxstyle='square', facecolor='white', alpha=1.)
ax.text(0.05, 0.95, r"$(a)$", transform=ax.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)


#plot 2
dum     = abs(phi_vor.flatten())
phi_max = max(dum)


axb  = axes([leftb, bottomb, widthb, heightb])
axb.set_xticks([])
axb.set_yticks([])
density = axb.pcolor(phi_vor, cmap='RdBu',vmin=-phi_max,vmax=phi_max)

props = dict(boxstyle='square', facecolor='white', alpha=1.)
axb.text(0.05, 0.95, r"$(b)$", transform=axb.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

#plot 3
dum     = abs(phi_turb.flatten())
phi_max = max(dum)


axc  = axes([leftc, bottomc, widthc, heightc])
axc.set_xticks([])
axc.set_yticks([])

density = axc.pcolor(phi_turb, cmap='RdBu',vmin=-phi_max,vmax=phi_max)
pos = axc.get_position().get_points()
cbx = fig.add_axes([pos[0,1]+.41, pos[0,1],  0.025,pos[1,1]-pos[0,1]])
cbar = colorbar(density, cax=cbx, orientation='vertical')
cbar.locator = matplotlib.ticker.FixedLocator([-phi_max/1.0,0,phi_max/1.0])
cbar.ax.yaxis.set_minor_locator(MultipleLocator(phi_max/2))
cbar.update_ticks()
cbar.ax.xaxis.set_ticks_position('default')
cbar.ax.set_yticklabels([r"$-1$",r"$0$",r"$1$"],ha='left',fontsize=15)
cbar.ax.set_xlabel(r"${\phi}$",fontsize=25, rotation=360)
cbar.ax.xaxis.set_label_coords(5.25,.66)

props = dict(boxstyle='square', facecolor='white', alpha=1.)
axc.text(0.05, 0.95, r"$(c)$", transform=axc.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)




#plot 4
vx,vy=vel_field(128,128,psi_lam)

dum     = abs(psi_turb.flatten())
psi_max = max(dum)
axd  = axes([leftd, bottomd, widthd, heightd])
axd.set_xticks([])
axd.set_yticks([])

vorticity = axd.pcolor(psi_lam, cmap='BrBG',vmin=-psi_max, vmax=psi_max)
stream_points = [[10,30],[5,30],[15,30],[20,30],[30,30],[50,30],[60,30],[70,30],[80,30],[90,60],[100,30],[120,30]]#np.array(zip(np.arange(0,62,.5), -np.arange(0,62,.5)))

cont=axd.streamplot(X,Y,vx,vy,density=0.5,color='black',linewidth=.5,arrowsize=.5,start_points=stream_points)
axd.text(219.5e-1, -21.95, r"$k^{'}=-0.003$", size=1.5*1e1)



#plot 5
vx,vy=vel_field(128,128,psi_vor)
dum     = abs(psi_turb.flatten())
psi_max = max(dum)
axe  = axes([lefte, bottome, widthe, heighte])
axe.set_xticks([])
axe.set_yticks([])

vorticity = axe.pcolor(psi_vor, cmap='BrBG',vmin=-psi_max, vmax=psi_max)
cont=axe.streamplot(X,Y,vx,vy,density=0.4,color='black',linewidth=.5,arrowsize=.5,broken_streamlines=False)
axe.text(219.5e-1, -21.95, r"$k^{'}=-0.0035$", size=1.5*1e1)



#plot 6
vx,vy=vel_field(128,128,psi_turb)
dum     = abs(psi_turb.flatten())
psi_max = max(dum)
axf  = axes([leftf, bottomf, widthf, heightf])
axf.set_xticks([])
axf.set_yticks([])

vorticity = axf.pcolor(psi_turb, cmap='BrBG',vmin=-psi_max, vmax=psi_max)
cont=axf.streamplot(X,Y,vx,vy,density=0.45,color='black',linewidth=.5,arrowsize=.5,broken_streamlines=False)
axf.text(219.5e-1, -21.95, r"$k^{'}=-0.06$", size=1.5*1e1)


pos = axf.get_position().get_points()
cbx = fig.add_axes([pos[0,1]+.875, pos[0,1],  0.025,pos[1,1]-pos[0,1]])
cbar = colorbar(vorticity, cax=cbx, orientation='vertical')
cbar.locator = matplotlib.ticker.FixedLocator([-psi_max/1.0,0,psi_max/1.0])
cbar.ax.xaxis.set_minor_locator(MultipleLocator(psi_max/2))
cbar.update_ticks()
cbar.ax.xaxis.set_ticks_position('bottom')
cbar.ax.set_yticklabels([r"$-4.7$",r"$0$",r"$4.7$"],ha='left',fontsize=15)
cbar.ax.set_xlabel(r"$\tilde{\psi} $",fontsize=25, rotation=360)
cbar.ax.xaxis.set_label_coords(5.25,.65)




fig.savefig('states_combined.png',bbox_inches='tight', dpi=500)
plt.show()

