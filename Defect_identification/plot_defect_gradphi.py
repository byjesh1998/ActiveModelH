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
with open('phi_46000000.txt','r') as f:
	listl0=[]
	for line in f:
		strip_lines=line.strip()
		listli=strip_lines.split()
		#print(listli)
		m=listl0.append(listli)
	#print(listl0)

#print(vx)		
a=np.asfarray(listl0,float)


with open('psi_46000000.txt','r') as f:
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
        
    return S11,S12, np.real(d2psi**2),rho, d2psi,vx,grad_phi,dxpsi,dypsi

s11,s12,grad,rho,v,vx,vy,fx,fy=energy_spectrum(nx,ny,a,b)

#S,R=np.meshgrid(s11,s12)
#print(s11.shape)
r=np.sqrt(s11**2 + s12**2)

lattice_angle= fns.smoothening(s12,s11,64,dx,2) #interchange
w=lattice_angle


dum     = abs(a.flatten())
epr_max = max(dum)



fig2 = plt.figure(2,figsize=(5,5))
ax2 = fig2.add_subplot(221)
epr=ax2.pcolor(a,cmap='seismic',vmin=-epr_max,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax2.spines["bottom"].set_linewidth(2)
ax2.spines["left"].set_linewidth(2)
ax2.spines["right"].set_linewidth(2)
ax2.spines["top"].set_linewidth(2)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_title(r"$+\frac{1}{2}$",fontsize=15)
ax2.set_xlim(6,21)   #1
ax2.set_ylim(94,109)


ax2.annotate("1", (7, 108), color='black', weight='bold', fontsize=12, ha='left', va='center')

ax3 = fig2.add_subplot(222)
epr=ax3.pcolor(a,cmap='seismic',vmin=-epr_max,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax3.spines["bottom"].set_linewidth(2)
ax3.spines["left"].set_linewidth(2)
ax3.spines["right"].set_linewidth(2)
ax3.spines["top"].set_linewidth(2)
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_title(r"$-\frac{1}{2}$",fontsize=15)
ax3.set_xlim(42,57) #2
ax3.set_ylim(97,112)

plt.colorbar(epr,label=r"$ \phi$")

ax3.annotate("2", (43, 111), color='black', weight='bold', fontsize=12, ha='left', va='center')



ax4 = fig2.add_subplot(223)
epr=ax4.pcolor(a,cmap='seismic',vmin=-epr_max,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax4.spines["bottom"].set_linewidth(2)
ax4.spines["left"].set_linewidth(2)
ax4.spines["right"].set_linewidth(2)
ax4.spines["top"].set_linewidth(2)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.set_xlim(112,127)#3
ax4.set_ylim(85,100)



ax4.annotate("3", (113, 99), color='black', weight='bold', fontsize=12, ha='left', va='center')


ax5 = fig2.add_subplot(224)
epr=ax5.pcolor(a,cmap='seismic',vmin=-epr_max,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax5.spines["bottom"].set_linewidth(2)
ax5.spines["left"].set_linewidth(2)
ax5.spines["right"].set_linewidth(2)
ax5.spines["top"].set_linewidth(2)
ax5.set_xticks([])
ax5.set_yticks([])
ax5.set_xlim(105,120) #4
ax5.set_ylim(25,40)

ax5.annotate("4", (106, 39), color='black', weight='bold', fontsize=12, ha='left', va='center')


plt.savefig('def-phi.png',dpi=500)



dum     = abs(vy.flatten())
epr_max = max(dum)

fig3 = plt.figure(3,figsize=(5,5))
ax2 = fig3.add_subplot(221)
epr=ax2.pcolor(vy,cmap='rainbow',vmin=0,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax2.spines["bottom"].set_linewidth(2)
ax2.spines["left"].set_linewidth(2)
ax2.spines["right"].set_linewidth(2)
ax2.spines["top"].set_linewidth(2)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_title(r"$+\frac{1}{2}$",fontsize=15)
ax2.set_xlim(6,21)   #1
ax2.set_ylim(94,109)

ax2.annotate("1", (7, 108), color='black', weight='bold', fontsize=12, ha='left', va='center')

ax3 = fig3.add_subplot(222)
epr=ax3.pcolor(vy,cmap='rainbow',vmin=0,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax3.spines["bottom"].set_linewidth(2)
ax3.spines["left"].set_linewidth(2)
ax3.spines["right"].set_linewidth(2)
ax3.spines["top"].set_linewidth(2)
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_title(r"$-\frac{1}{2}$",fontsize=15)
ax3.set_xlim(42,57) #2
ax3.set_ylim(97,112)

plt.colorbar(epr,label=r"$\nabla \phi$")

ax3.annotate("2", (43, 111), color='black', weight='bold', fontsize=12, ha='left', va='center')



ax4 = fig3.add_subplot(223)
epr=ax4.pcolor(vy,cmap='rainbow',vmin=0,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax4.spines["bottom"].set_linewidth(2)
ax4.spines["left"].set_linewidth(2)
ax4.spines["right"].set_linewidth(2)
ax4.spines["top"].set_linewidth(2)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.set_xlim(112,127)#3
ax4.set_ylim(85,100)

ax4.annotate("3", (113, 99), color='black', weight='bold', fontsize=12, ha='left', va='center')


ax5 = fig3.add_subplot(224)
epr=ax5.pcolor(vy,cmap='rainbow',vmin=0,vmax=epr_max)
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax5.spines["bottom"].set_linewidth(2)
ax5.spines["left"].set_linewidth(2)
ax5.spines["right"].set_linewidth(2)
ax5.spines["top"].set_linewidth(2)
ax5.set_xticks([])
ax5.set_yticks([])
ax5.set_xlim(105,120) #4
ax5.set_ylim(25,40)

ax5.annotate("4", (106, 39), color='black', weight='bold', fontsize=12, ha='left', va='center')


plt.savefig('def_grad_phi.png',dpi=500)

import matplotlib.colors as colors


fig4 = plt.figure(4,figsize=(5,5))
ax2 = fig4.add_subplot(221)
epr=ax2.pcolor(vy,cmap='rainbow',norm=colors.LogNorm(vmin=1e-4,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax2.spines["bottom"].set_linewidth(2)
ax2.spines["left"].set_linewidth(2)
ax2.spines["right"].set_linewidth(2)
ax2.spines["top"].set_linewidth(2)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_title(r"$+\frac{1}{2}$",fontsize=15)
ax2.set_xlim(6,21)   #1
ax2.set_ylim(94,109)

ax2.annotate("1", (7, 108), color='black', weight='bold', fontsize=12, ha='left', va='center')

ax3 = fig4.add_subplot(222)
epr=ax3.pcolor(vy,cmap='rainbow',norm=colors.LogNorm(vmin=1e-4,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax3.spines["bottom"].set_linewidth(2)
ax3.spines["left"].set_linewidth(2)
ax3.spines["right"].set_linewidth(2)
ax3.spines["top"].set_linewidth(2)
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_title(r"$-\frac{1}{2}$",fontsize=15)
ax3.set_xlim(42,57) #2
ax3.set_ylim(97,112)

plt.colorbar(epr,label=r"$\nabla \phi$")
plt.title("12")

ax3.annotate("2", (43, 111), color='black', weight='bold', fontsize=12, ha='left', va='center')



ax4 = fig4.add_subplot(223)
epr=ax4.pcolor(vy,cmap='rainbow',norm=colors.LogNorm(vmin=1e-4,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax4.spines["bottom"].set_linewidth(2)
ax4.spines["left"].set_linewidth(2)
ax4.spines["right"].set_linewidth(2)
ax4.spines["top"].set_linewidth(2)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.set_xlim(112,127)#3
ax4.set_ylim(85,100)

ax4.annotate("3", (113, 99), color='black', weight='bold', fontsize=12, ha='left', va='center')


ax5 = fig4.add_subplot(224)
epr=ax5.pcolor(vy,cmap='rainbow',norm=colors.LogNorm(vmin=1e-4,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
ax5.spines["bottom"].set_linewidth(2)
ax5.spines["left"].set_linewidth(2)
ax5.spines["right"].set_linewidth(2)
ax5.spines["top"].set_linewidth(2)
ax5.set_xticks([])
ax5.set_yticks([])
ax5.set_xlim(105,120) #4
ax5.set_ylim(25,40)

ax5.annotate("4", (106, 39), color='black', weight='bold', fontsize=12, ha='left', va='center')


plt.savefig('def_grad_phi2.png',dpi=500)


fig5 = plt.figure(5,figsize=(5,5))
ax2 = fig5.add_subplot(221)
#lattice_angle= fns.smoothening(fx,fy,64,dx,2)*2 #interchange
#l=lattice_angle

#epr=ax2.pcolor(vy,cmap='rainbow',norm=colors.LogNorm(vmin=1e-4,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
plt.quiver(X,Y,fx,fy,scale=8., color='red',width=0.01,headwidth=0.8, headlength=0.8)
ax2.spines["bottom"].set_linewidth(2)
ax2.spines["left"].set_linewidth(2)
ax2.spines["right"].set_linewidth(2)
ax2.spines["top"].set_linewidth(2)
ax2.set_xticks([])
ax2.set_yticks([])
ax2.set_title("nematic director fields",fontsize=15)
ax2.set_xlim(6,21)   #1
ax2.set_ylim(94,109)

ax2.annotate("1", (7, 108), color='black', weight='bold', fontsize=12, ha='left', va='center')

ax3 = fig5.add_subplot(222)
#epr=ax3.pcolor(vy,cmap='rainbow',norm=colors.LogNorm(vmin=1e-4,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
         
plt.quiver(X,Y,fx,fy,scale=8., color='red',width=0.01,headwidth=0.8, headlength=0.8)
ax3.spines["bottom"].set_linewidth(2)
ax3.spines["left"].set_linewidth(2)
ax3.spines["right"].set_linewidth(2)
ax3.spines["top"].set_linewidth(2)
ax3.set_xticks([])
ax3.set_yticks([])
ax3.set_title(r"$-\frac{1}{2}$",fontsize=15)
ax3.set_xlim(42,57) #2
ax3.set_ylim(97,112)

#plt.colorbar(epr,label=r"$\nabla \phi$")
plt.title(r"$\nabla \phi$",color='red')

ax3.annotate("2", (43, 111), color='black', weight='bold', fontsize=12, ha='left', va='center')



ax4 = fig5.add_subplot(223)
#epr=ax4.pcolor(vy,cmap='rainbow',norm=colors.LogNorm(vmin=1e-4,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
plt.quiver(X,Y,fx,fy,scale=8., color='red',width=0.01,headwidth=0.8, headlength=0.8)
ax4.spines["bottom"].set_linewidth(2)
ax4.spines["left"].set_linewidth(2)
ax4.spines["right"].set_linewidth(2)
ax4.spines["top"].set_linewidth(2)
ax4.set_xticks([])
ax4.set_yticks([])
ax4.set_xlim(112,127)#3
ax4.set_ylim(85,100)

ax4.annotate("3", (113, 99), color='black', weight='bold', fontsize=12, ha='left', va='center')


ax5 = fig5.add_subplot(224)
#epr=ax5.pcolor(vy,cmap='rainbow',norm=colors.LogNorm(vmin=1e-4,vmax=epr_max))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=18., color='black',
         width=0.01, headwidth=0.8, headlength=0.8)
         
         
plt.quiver(X,Y,fx,fy,scale=8., color='red',width=0.01,headwidth=0.8, headlength=0.8)
ax5.spines["bottom"].set_linewidth(2)
ax5.spines["left"].set_linewidth(2)
ax5.spines["right"].set_linewidth(2)
ax5.spines["top"].set_linewidth(2)
ax5.set_xticks([])
ax5.set_yticks([])
ax5.set_xlim(105,120) #4
ax5.set_ylim(25,40)

ax5.annotate("4", (106, 39), color='black', weight='bold', fontsize=12, ha='left', va='center')


plt.savefig('def_grad_phi_dir.png',dpi=500)



fig6 = plt.figure(6,figsize=(5,5))


from mpl_toolkits import mplot3d


ax = plt.axes(projection='3d')
ax.plot_surface(X,Y,a,cmap='viridis')

plt.savefig('phi_3d.png',dpi=500)
plt.show()





    
