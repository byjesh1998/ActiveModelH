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
        
    return GK,HK,S11,S12

GK,HK,s11,s12=energy_spectrum(nx,ny,a,b)

#S,R=np.meshgrid(s11,s12)
#print(s11.shape)
#r=np.sqrt(s11**2 + s12**2)

lattice_angle= fns.smoothening(s12,s11,64,dx,4) #interchange
w=lattice_angle



###########################################################################
#first plot
#GK=HK
plt.pcolor(GK,cmap='seismic')
#plt.contour(X,Y,GK,20)
plt.colorbar()

plt.quiver(X,Y,np.cos(w),np.sin(w),scale=10)


plt.show()






    
