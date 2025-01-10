import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import matplotlib.animation as animation
#plt.rcParams['animation.ffmpeg_path'] = r"C:\some_path\ffmpeg.exe"   # if necessary
import glob, os
import re
#import pyfftw
import pandas as pd
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Rectangle
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as mticker
from matplotlib import rc
import fns
import matplotlib as mpl
import matplotlib.colors as colors

import fns

f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.1e' % x))
fmt = mticker.FuncFormatter(g)



nx=int(128)
ny=int(128)
Nx=nx
x=np.linspace(0,nx,nx)
y=np.linspace(0,nx,nx)

X,Y=np.meshgrid(x,y)

#vx = pyfftw.empty_aligned((nx,ny),dtype= 'complex128')
#vy=np.zeros((nx,ny))
#a=[]
#b=[]
with open('phi_turb.txt','r') as f:
	listl0=[]
	for line in f:
		strip_lines=line.strip()
		listli=strip_lines.split()
		#print(listli)
		m=listl0.append(listli)
	#print(listl0)

#print(vx)		
a=np.asfarray(listl0,float)


with open('psi_turb.txt','r') as f:
	list10=[]
	for line in f:
		strip_lines=line.strip()
		listli=strip_lines.split()
		#print(listli)
		m=list10.append(listli)
	#print(listl0)

#print(vx)		
b=np.asfarray(list10,float)

c=a
#vy=b


dx=0.5000
# compute the energy spectrum numerically
def energy_spectrum(nx,ny,psi,b):
    
    '''
    Computation of energy spectrum and maximum wavenumber from vorticity field
    
    Inputs
    ------
    nx,ny : number of grid points in x and y direction
    w : vorticity field in physical spce (including periodic boundaries)
    
    Output
    ------
    en : energy spectrum computed from vorticity field
    n : maximum wavenumber
    '''
    
    epsilon = 1.0e-8

    kx = np.empty(nx)
    ky = np.empty(ny)
    
    kx[0:int(nx/2)] = np.float64(np.arange(0,int(nx/2)))* 2*np.pi/(np.float64(nx)*dx)
    kx[int(nx/2):nx] = np.float64(np.arange(-int(nx/2),0))* 2*np.pi/(np.float64(nx)*dx)

    ky[0:ny] = kx[0:ny]
    
    kx[0] = epsilon
    ky[0] = epsilon
    vx=np.random.randint(0,10,(nx,ny))
    ky, kx = np.meshgrid(kx, ky, indexing='ij')
    
    #a = pyfftw.empty_aligned((nx,ny),dtype= 'complex128')
    #b = pyfftw.empty_aligned((nx,ny),dtype= 'complex128')

    #fft_object = pyfftw.FFTW(a, b, axes = (0,1), direction = 'FFTW_FORWARD')
    
    #wf = fft_object(w[0:nx,0:ny]) 
    #vxf = fft_object(vx[0:nx,0:ny] )
    #vyf = fft_object(vy[0:nx,0:ny]) 
    
    psif=np.fft.fft2(psi)
    
    #print(len(vxf))
    
    
    kk = np.sqrt(kx[:,:]**2 + ky[:,:]**2)
    kx2=kx[:,:]**2
    ky2=ky[:,:]**2
    #s11=0.5*(- ky2*psif + kx2 * psif)
    #s12=  -(1j*kx[:,:] * psif) * (1j*ky[:,:]* psif)
    dypsi_ft=(1j*ky[:,:]* psif)
    dxpsi_ft=(1j*kx[:,:]* psif)
    
 
    
    dxpsi=np.real(np.fft.ifft2(dxpsi_ft))

    dypsi=np.real(np.fft.ifft2(dypsi_ft))

    dyypsi_ft=(-ky[:,:]**2* psif)
    dxxpsi_ft=(-kx[:,:]**2* psif)

    dxypsi_ft=(-kx[:,:]*ky[:,:]* psif)
       
    dxxpsi=np.real(np.fft.ifft2(dxxpsi_ft))
    dyypsi=np.real(np.fft.ifft2(dyypsi_ft))
    dxypsi=np.real(np.fft.ifft2(dxypsi_ft))

    GK=(dxxpsi*dyypsi-dxypsi**2)/((1+dxpsi**2+dypsi**2)**2)

    HK=((1+dxpsi**2)*dyypsi +(1+dypsi**2)*dxxpsi - 2*dxpsi*dypsi*dxypsi)/(2*(dxpsi**2+dypsi**2+1)**(1.5))

    grad_phi=np.sqrt(dxpsi **2+dypsi **2)


    S11=0.5 * ((dypsi **2)-(dxpsi)**2) #add -
    S12= - dxpsi * dypsi
    
    r=S11**2+S12**2
    
    s11_ft=np.fft.fft2(S11)
    s12_ft=np.fft.fft2(S12)
    
    dxs11_ft=(1j*kx[:,:]* s11_ft)
    dys11_ft=(1j*ky[:,:]* s11_ft)
    
    dxs12_ft=(1j*kx[:,:]* s12_ft)
    dys12_ft=(1j*ky[:,:]* s12_ft)
    
    dxs11=np.real(np.fft.ifft2(dxs11_ft))
    dys11=np.real(np.fft.ifft2(dys11_ft))
    dxs12=np.real(np.fft.ifft2(dxs12_ft))
    dys12=np.real(np.fft.ifft2(dys12_ft))
    
    rho=dxs11*dys12-dxs12*dys11
    
    
    bf=np.fft.fft2(b)
    d2psi_ft=  -kk**2 *bf
    d4psi_ft=  kk**4 *bf
    vx_ft= 1j*ky[:,:]* bf
    vy_ft= -1j*kx[:,:]* bf


    vx=np.real(np.fft.ifft2(vx_ft))
    vy=np.real(np.fft.ifft2(vy_ft))
    
    
    d2psi=np.real(np.fft.ifft2(d2psi_ft))
    d4psi=np.real(np.fft.ifft2(d4psi_ft))
    
    #grad term
    
    d2dxpsi_ft= - kk **2 *dxpsi_ft
    d2dypsi_ft= - kk **2 *dypsi_ft
    
    d2dxpsi=np.fft.ifft2(d2dxpsi_ft)
    d2dypsi=np.fft.ifft2(d2dypsi_ft)
    
    grad=(dypsi * d2dxpsi - dxpsi * d2dypsi) *b#*(0.2)
        
    return S11,S12, np.real(d2psi**2),rho, d2psi,vx,grad_phi,dxpsi,dypsi,GK,HK

s11,s12,grad,rho,v,vx,vy,fx,fy,GK,HK=energy_spectrum(nx,ny,a,b)

#S,R=np.meshgrid(s11,s12)
#print(s11.shape)
r=np.sqrt(s11**2 + s12**2)

lattice_angle= fns.smoothening(s12,s11,64,dx,4) #interchange
w=lattice_angle



###########################################################################
#first plot
grad=vy
dum     = abs(grad.flatten())
epr_max = max(dum)
	

params = {#'legend.fontsize': 6,
          #'text.latex.preamble': [r"\usepackage{amstext}",],
          #'font',**{'family':'sans-serif','sans-serif':['Helvetica']},
          'axes.linewidth': 2.5,
          #'axes.labelsize': 6.5,
          #'text.fontsize': 4,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,	
          'text.usetex': True
          #'figure.figsize': fig_size
          }
mpl.rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

###########################################################################


#second plot
	
	
fig2 = plt.figure(2,figsize=(5,5))
ax2 = fig2.add_subplot(221)
epr=ax2.pcolor(grad,cmap='rainbow',vmin=0,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax2.spines["bottom"].set_linewidth(2)
ax2.spines["left"].set_linewidth(2)
ax2.spines["right"].set_linewidth(2)
ax2.spines["top"].set_linewidth(2)
ax2.set_xticks([])
ax2.set_yticks([])
#ax2.set_title(r"$+\frac{1}{2}$",fontsize=15)
ax2.set_xlim(18,32)   #1
ax2.set_ylim(66,81)

ax2.annotate("1", (19, 79), color='black', weight='bold', fontsize=22, ha='left', va='center')

ax3 = fig2.add_subplot(222)
epr=ax3.pcolor(grad,cmap='rainbow',vmin=0,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax3.spines["bottom"].set_linewidth(2)
ax3.spines["left"].set_linewidth(2)
ax3.spines["right"].set_linewidth(2)
ax3.spines["top"].set_linewidth(2)
ax3.set_xticks([])
ax3.set_yticks([])
#ax3.set_title(r"$-\frac{1}{2}$",fontsize=15)
ax3.set_xlim(39,54) #2
ax3.set_ylim(69,84)

ax3.annotate("2", (40, 82), color='black', weight='bold', fontsize=22, ha='left', va='center')



ax4 = fig2.add_subplot(223)
epr=ax4.pcolor(grad,cmap='rainbow',vmin=0,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax4.spines["bottom"].set_linewidth(2)
ax4.spines["left"].set_linewidth(2)
ax4.spines["right"].set_linewidth(2)
ax4.spines["top"].set_linewidth(2)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.set_xlim(82,97)#3
ax4.set_ylim(78,93)

ax4.annotate("3", (83, 91), color='black', weight='bold', fontsize=22, ha='left', va='center')
ax4.set_xlabel(r"$+\frac{1}{2}$ $defect$",fontsize=25)

ax5 = fig2.add_subplot(224)
epr=ax5.pcolor(grad,cmap='rainbow',vmin=0,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax5.spines["bottom"].set_linewidth(2)
ax5.spines["left"].set_linewidth(2)
ax5.spines["right"].set_linewidth(2)
ax5.spines["top"].set_linewidth(2)
ax5.set_xticks([])
ax5.set_yticks([])
ax5.set_xlim(72,87) #4
ax5.set_ylim(5,20)

ax5.annotate("4", (73, 18), color='black', weight='bold', fontsize=22, ha='left', va='center')

ax5.set_xlabel(r"$-\frac{1}{2}$ $defect$",fontsize=25)
pos = ax2.get_position().get_points()
pos1 = ax3.get_position().get_points()
cax = fig2.add_axes([
    pos[0,0], pos[1,1]+0.02, pos1[1,0]-pos[0,0], 0.025
])
cbar = fig2.colorbar(epr, cax=cax, orientation='horizontal')
cbar.locator = matplotlib.ticker.FixedLocator([0,epr_max/4,epr_max/2,3*epr_max/4,epr_max])
cbar.ax.xaxis.set_minor_locator(MultipleLocator(1/2))
cbar.update_ticks()
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.set_xticklabels([r"$0$",r"$0.2$",r"$0.4$",r"$0.6$",r"$0.8$"],fontsize=22)
#cbar.ax.set_xticklabels([r"$-3 \times 10^{-5}$", r"$0$",r"$3 \times 10^{-5}$"])
#cbar.ax.set_xticklabels(["-{}".format(fmt(epr_max)),'0',"{}".format(fmt(epr_max))],ha='left')
cbar.ax.set_ylabel(r"${\nabla \phi}$",fontsize=30, rotation=360)
cbar.ax.yaxis.set_label_coords(-0.08,-1.5)

plt.subplots_adjust(wspace=0.051, hspace=0.051)


plt.savefig('def_gradphi.png',dpi=500)


###########################################################################
	

#ax2.set_aspect('equal')

###########################################################################


#third plot
dum     = abs(a.flatten())
v_max = max(dum)	
	
fig3 = plt.figure(3,figsize=(5,5))
ax = fig3.add_subplot(111)
vorticity=ax.pcolor(a,cmap='RdBu',vmin=-v_max,vmax=v_max)
#plt.quiver(X,Y,np.cos(w),np.sin(w),scale=20., color='black',
#         width=0.005, headwidth=0.8, headlength=0.8)
#plt.streamplot(X,Y,-vy,-vx, color='white',arrowsize=1,arrowstyle='-|>',density=1)
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["top"].set_linewidth(2)
ax.set_xticks([])
ax.set_yticks([])

pos = ax.get_position().get_points()
cax = fig3.add_axes([
    pos[0,0], pos[1,1]+0.01, pos[1,0]-pos[0,0], 0.025
])
cbar = fig3.colorbar(vorticity, cax=cax, orientation='horizontal')
cbar.locator = matplotlib.ticker.FixedLocator([-1,0,1])
cbar.ax.xaxis.set_minor_locator(MultipleLocator(1/2))
cbar.update_ticks()
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.set_xticklabels([r"$-1$",r"$0$",r"$1$"],fontsize=22)
#cbar.ax.set_xticklabels([r"$-3 \times 10^{-5}$", r"$0$",r"$3 \times 10^{-5}$"])
#cbar.ax.set_xticklabels(["-{}".format(fmt(epr_max)),'0',"{}".format(fmt(epr_max))],ha='left')
cbar.ax.set_ylabel(r"${\phi}$",fontsize=30, rotation=360)
cbar.ax.yaxis.set_label_coords(-0.08,-1.5)

ax.add_patch(Rectangle((39, 69), 15, 15,edgecolor='black',
                    facecolor='none',
                    lw=2.5))

ax.annotate("2", (39, 89), color='black', weight='bold', fontsize=25, ha='left', va='center')



ax.add_patch(Rectangle((18, 66), 15, 15,edgecolor='black',
                    facecolor='none',
                    lw=2.5))
ax.annotate("1", (18, 86), color='black', weight='bold', fontsize=25, ha='left', va='center')



ax.add_patch(Rectangle((82, 78), 15, 15,edgecolor='black',
                    facecolor='none',
                    lw=2.5))
ax.annotate("3", (82, 98), color='black', weight='bold', fontsize=25, ha='left', va='center')


ax.add_patch(Rectangle((72, 5), 15, 15,edgecolor='black',
                    facecolor='none',
                    lw=2.5))
ax.annotate("4", (72, 25), color='black', weight='bold', fontsize=25, ha='left', va='center')


plt.savefig('phi.png',dpi=500)

	



##################################################################################3

fig4 = plt.figure(4,figsize=(5,5))
ax2 = fig4.add_subplot(221)
epr=ax2.pcolor(grad,cmap='rainbow',norm=colors.LogNorm(vmin=1e-3,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax2.spines["bottom"].set_linewidth(2)
ax2.spines["left"].set_linewidth(2)
ax2.spines["right"].set_linewidth(2)
ax2.spines["top"].set_linewidth(2)
ax2.set_xticks([])
ax2.set_yticks([])
#ax2.set_title(r"$+\frac{1}{2}$",fontsize=15)
ax2.set_xlim(18,32)   #1
ax2.set_ylim(66,81)

ax2.annotate("1", (19, 79), color='black', weight='bold', fontsize=22, ha='left', va='center')

ax3 = fig4.add_subplot(222)
epr=ax3.pcolor(grad,cmap='rainbow',norm=colors.LogNorm(vmin=1e-3,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax3.spines["bottom"].set_linewidth(2)
ax3.spines["left"].set_linewidth(2)
ax3.spines["right"].set_linewidth(2)
ax3.spines["top"].set_linewidth(2)
ax3.set_xticks([])
ax3.set_yticks([])
#ax3.set_title(r"$-\frac{1}{2}$",fontsize=15)
ax3.set_xlim(39,54) #2
ax3.set_ylim(69,84)

ax3.annotate("2", (40, 82), color='black', weight='bold', fontsize=22, ha='left', va='center')



ax4 = fig4.add_subplot(223)
epr=ax4.pcolor(grad,cmap='rainbow',norm=colors.LogNorm(vmin=1e-3,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax4.spines["bottom"].set_linewidth(2)
ax4.spines["left"].set_linewidth(2)
ax4.spines["right"].set_linewidth(2)
ax4.spines["top"].set_linewidth(2)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.set_xlim(82,97)#3
ax4.set_ylim(78,93)

ax4.annotate("3", (83, 91), color='black', weight='bold', fontsize=22, ha='left', va='center')
ax4.set_xlabel(r"$+\frac{1}{2}$ $defect$",fontsize=25)

ax5 = fig4.add_subplot(224)
epr=ax5.pcolor(grad,cmap='rainbow',norm=colors.LogNorm(vmin=1e-3,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax5.spines["bottom"].set_linewidth(2)
ax5.spines["left"].set_linewidth(2)
ax5.spines["right"].set_linewidth(2)
ax5.spines["top"].set_linewidth(2)
ax5.set_xticks([])
ax5.set_yticks([])
ax5.set_xlim(72,87) #4
ax5.set_ylim(5,20)

ax5.annotate("4", (73, 18), color='black', weight='bold', fontsize=22, ha='left', va='center')

ax5.set_xlabel(r"$-\frac{1}{2}$ $defect$",fontsize=25)
pos = ax2.get_position().get_points()
pos1 = ax3.get_position().get_points()
cax = fig4.add_axes([
    pos[0,0], pos[1,1]+0.02, pos1[1,0]-pos[0,0], 0.025
])
cbar = fig4.colorbar(epr, cax=cax, orientation='horizontal')
##cbar.locator = matplotlib.ticker.FixedLocator([0,epr_max/4,epr_max/2,3*epr_max/4,epr_max])
##cbar.ax.xaxis.set_minor_locator(MultipleLocator(1/2))
cbar.update_ticks()
cbar.ax.xaxis.set_ticks_position('top')
##cbar.ax.set_xticklabels([r"$0$",r"$0.2$",r"$0.4$",r"$0.6$",r"$0.8$"],fontsize=22)
#cbar.ax.set_xticklabels([r"$-3 \times 10^{-5}$", r"$0$",r"$3 \times 10^{-5}$"])
#cbar.ax.set_xticklabels(["-{}".format(fmt(epr_max)),'0',"{}".format(fmt(epr_max))],ha='left')
cbar.ax.set_ylabel(r"${\nabla \phi}$",fontsize=30, rotation=360)
cbar.ax.yaxis.set_label_coords(-0.08,-1.5)

plt.subplots_adjust(wspace=0.051, hspace=0.051)

plt.savefig('def_gradphi_2.png',dpi=500)





#plt.figure(5)
fig5 = plt.figure(1,figsize=(5,5))
ax5 = fig5.add_subplot(111)
epr=ax5.pcolor(grad,cmap='rainbow',norm=colors.LogNorm(vmin=1e-4,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax5.spines["bottom"].set_linewidth(2)
ax5.spines["left"].set_linewidth(2)
ax5.spines["right"].set_linewidth(2)
ax5.spines["top"].set_linewidth(2)
ax5.set_xticks([])
ax5.set_yticks([])

pos = ax5.get_position().get_points()

cax = fig5.add_axes([
    pos[0,0], pos[1,1]+0.02, pos[1,0]-pos[0,0], 0.025
])
cbar = fig5.colorbar(epr, cax=cax, orientation='horizontal')
##cbar.locator = matplotlib.ticker.FixedLocator([0,epr_max/4,epr_max/2,3*epr_max/4,epr_max])
##cbar.ax.xaxis.set_minor_locator(MultipleLocator(1/2))
cbar.update_ticks()
cbar.ax.xaxis.set_ticks_position('top')
##cbar.ax.set_xticklabels([r"$0$",r"$0.2$",r"$0.4$",r"$0.6$",r"$0.8$"],fontsize=22)
#cbar.ax.set_xticklabels([r"$-3 \times 10^{-5}$", r"$0$",r"$3 \times 10^{-5}$"])
#cbar.ax.set_xticklabels(["-{}".format(fmt(epr_max)),'0',"{}".format(fmt(epr_max))],ha='left')
cbar.ax.set_ylabel(r"${\nabla \phi}$",fontsize=30, rotation=360)
cbar.ax.yaxis.set_label_coords(-0.08,-1.5)
#ax5.set_xlim(72,87) #4
#ax5.set_ylim(5,20)

###########################################################################


#second plot
	
grad=HK
dum     = abs(grad.flatten())
epr_max = max(dum)	

fig6 = plt.figure(6,figsize=(5,5))
ax2 = fig6.add_subplot(221)
epr=ax2.pcolor(grad,cmap='bwr',vmin=-epr_max,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax2.spines["bottom"].set_linewidth(2)
ax2.spines["left"].set_linewidth(2)
ax2.spines["right"].set_linewidth(2)
ax2.spines["top"].set_linewidth(2)
ax2.set_xticks([])
ax2.set_yticks([])
#ax2.set_title(r"$+\frac{1}{2}$",fontsize=15)
ax2.set_xlim(18,32)   #1
ax2.set_ylim(66,81)

ax2.annotate("1", (19, 79), color='black', weight='bold', fontsize=22, ha='left', va='center')

ax3 = fig6.add_subplot(222)
epr=ax3.pcolor(grad,cmap='bwr',vmin=-epr_max,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax3.spines["bottom"].set_linewidth(2)
ax3.spines["left"].set_linewidth(2)
ax3.spines["right"].set_linewidth(2)
ax3.spines["top"].set_linewidth(2)
ax3.set_xticks([])
ax3.set_yticks([])
#ax3.set_title(r"$-\frac{1}{2}$",fontsize=15)
ax3.set_xlim(39,54) #2
ax3.set_ylim(69,84)

ax3.annotate("2", (40, 82), color='black', weight='bold', fontsize=22, ha='left', va='center')



ax4 = fig6.add_subplot(223)
epr=ax4.pcolor(grad,cmap='bwr',vmin=-epr_max,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax4.spines["bottom"].set_linewidth(2)
ax4.spines["left"].set_linewidth(2)
ax4.spines["right"].set_linewidth(2)
ax4.spines["top"].set_linewidth(2)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.set_xlim(82,97)#3
ax4.set_ylim(78,93)

ax4.annotate("3", (83, 91), color='black', weight='bold', fontsize=22, ha='left', va='center')
ax4.set_xlabel(r"$+\frac{1}{2}$ $defect$",fontsize=25)

ax5 = fig6.add_subplot(224)
epr=ax5.pcolor(grad,cmap='bwr',vmin=-epr_max,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax5.spines["bottom"].set_linewidth(2)
ax5.spines["left"].set_linewidth(2)
ax5.spines["right"].set_linewidth(2)
ax5.spines["top"].set_linewidth(2)
ax5.set_xticks([])
ax5.set_yticks([])
ax5.set_xlim(72,87) #4
ax5.set_ylim(5,20)

ax5.annotate("4", (73, 18), color='black', weight='bold', fontsize=22, ha='left', va='center')

ax5.set_xlabel(r"$-\frac{1}{2}$ $defect$",fontsize=25)
pos = ax2.get_position().get_points()
pos1 = ax3.get_position().get_points()
cax = fig6.add_axes([
    pos[0,0], pos[1,1]+0.02, pos1[1,0]-pos[0,0], 0.025
])
cbar = fig6.colorbar(epr, cax=cax, orientation='horizontal')
cbar.locator = matplotlib.ticker.FixedLocator([-epr_max,-epr_max/2,0,epr_max/2,epr_max])
cbar.ax.xaxis.set_minor_locator(MultipleLocator(1/2))
cbar.update_ticks()
cbar.ax.xaxis.set_ticks_position('top')
cbar.ax.set_xticklabels([r"$-0.25$",r"$-0.125$",r"$0$",r"$0.125$",r"$0.25$"],fontsize=22)
#cbar.ax.set_xticklabels([r"$-3 \times 10^{-5}$", r"$0$",r"$3 \times 10^{-5}$"])
#cbar.ax.set_xticklabels(["-{}".format(fmt(epr_max)),'0',"{}".format(fmt(epr_max))],ha='left')
cbar.ax.set_ylabel(r"${H}$",fontsize=30, rotation=360)
cbar.ax.yaxis.set_label_coords(-0.08,-1.5)

plt.subplots_adjust(wspace=0.051, hspace=0.051)

plt.savefig('def_mean_curv.png',dpi=500)



###########################################################################


#second plot
	
grad=HK
dum     = abs(grad.flatten())
epr_max = max(dum)	

fig7 = plt.figure(7,figsize=(5,5))
ax2 = fig7.add_subplot(221)
#epr=ax2.pcolor(grad,cmap='bwr',vmin=-epr_max,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
plt.quiver(X,Y,fx,fy,scale=6., color='red',width=0.01)
ax2.spines["bottom"].set_linewidth(2)
ax2.spines["left"].set_linewidth(2)
ax2.spines["right"].set_linewidth(2)
ax2.spines["top"].set_linewidth(2)
ax2.set_xticks([])
ax2.set_yticks([])
#ax2.set_title(r"$+\frac{1}{2}$",fontsize=15)
ax2.set_xlim(18,32)   #1
ax2.set_ylim(66,81)

ax2.annotate("1", (19, 79), color='black', weight='bold', fontsize=22, ha='left', va='center')

ax3 = fig7.add_subplot(222)
#epr=ax3.pcolor(grad,cmap='bwr',vmin=-epr_max,vmax=epr_max)
plt.quiver(X,Y,fx,fy,scale=6., color='red',width=0.01)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax3.spines["bottom"].set_linewidth(2)
ax3.spines["left"].set_linewidth(2)
ax3.spines["right"].set_linewidth(2)
ax3.spines["top"].set_linewidth(2)
ax3.set_xticks([])
ax3.set_yticks([])
#ax3.set_title(r"$-\frac{1}{2}$",fontsize=15)
ax3.set_xlim(39,54) #2
ax3.set_ylim(69,84)

ax3.annotate("2", (40, 82), color='black', weight='bold', fontsize=22, ha='left', va='center')



ax4 = fig7.add_subplot(223)
#epr=ax4.pcolor(grad,cmap='bwr',vmin=-epr_max,vmax=epr_max)
plt.quiver(X,Y,fx,fy,scale=6., color='red',width=0.01)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax4.spines["bottom"].set_linewidth(2)
ax4.spines["left"].set_linewidth(2)
ax4.spines["right"].set_linewidth(2)
ax4.spines["top"].set_linewidth(2)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.set_xlim(82,97)#3
ax4.set_ylim(78,93)

ax4.annotate("3", (83, 91), color='black', weight='bold', fontsize=22, ha='left', va='center')
ax4.set_xlabel(r"$+\frac{1}{2}$ $defect$",fontsize=25)

ax5 = fig7.add_subplot(224)
#epr=ax5.pcolor(grad,cmap='bwr',vmin=-epr_max,vmax=epr_max)
plt.quiver(X,Y,fx,fy,scale=6., color='red',width=0.01)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=14., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax5.spines["bottom"].set_linewidth(2)
ax5.spines["left"].set_linewidth(2)
ax5.spines["right"].set_linewidth(2)
ax5.spines["top"].set_linewidth(2)
ax5.set_xticks([])
ax5.set_yticks([])
ax5.set_xlim(72,87) #4
ax5.set_ylim(5,20)

ax5.annotate("4", (73, 18), color='black', weight='bold', fontsize=22, ha='left', va='center')

ax5.set_xlabel(r"$-\frac{1}{2}$ $defect$",fontsize=25)
pos = ax2.get_position().get_points()
pos1 = ax3.get_position().get_points()

plt.subplots_adjust(wspace=0.051, hspace=0.051)
# cax = fig6.add_axes([
#     pos[0,0], pos[1,1]+0.02, pos1[1,0]-pos[0,0], 0.025
# ])
# cbar = fig6.colorbar(epr, cax=cax, orientation='horizontal')
# cbar.locator = matplotlib.ticker.FixedLocator([-epr_max,-epr_max/2,0,epr_max/2,epr_max])
# cbar.ax.xaxis.set_minor_locator(MultipleLocator(1/2))
# cbar.update_ticks()
# cbar.ax.xaxis.set_ticks_position('top')
# cbar.ax.set_xticklabels([r"$-0.25$",r"$-0.125$",r"$0$",r"$0.125$",r"$0.25$"],fontsize=22)
# #cbar.ax.set_xticklabels([r"$-3 \times 10^{-5}$", r"$0$",r"$3 \times 10^{-5}$"])
# #cbar.ax.set_xticklabels(["-{}".format(fmt(epr_max)),'0',"{}".format(fmt(epr_max))],ha='left')
# cbar.ax.set_ylabel(r"${H}$",fontsize=30, rotation=360)
# cbar.ax.yaxis.set_label_coords(-0.08,-1.5)

plt.savefig('def_fx.png',dpi=500)


plt.show()






    
