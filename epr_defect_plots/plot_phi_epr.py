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
    
    
    d2psi=np.real(np.fft.ifft2(d2psi_ft))/((0.06/1.67))
    d4psi=np.real(np.fft.ifft2(d4psi_ft))
    
    #grad term
    
    d2dxpsi_ft= - kk **2 *dxpsi_ft
    d2dypsi_ft= - kk **2 *dypsi_ft
    
    d2dxpsi=np.fft.ifft2(d2dxpsi_ft)
    d2dypsi=np.fft.ifft2(d2dypsi_ft)
    
    grad=d2psi **2*1.67#(dypsi * d2dxpsi - dxpsi * d2dypsi) *b#*(0.2)
        
    return S11,S12, np.real(d2psi**2),rho, d2psi,vx,grad_phi,dxpsi,dypsi

s11,s12,grad,rho,v,vx,vy,fx,fy=energy_spectrum(nx,ny,a,b)

#S,R=np.meshgrid(s11,s12)
#print(s11.shape)
r=np.sqrt(s11**2 + s12**2)

lattice_angle= fns.smoothening(s12,s11,64,dx,6) #interchange
w=lattice_angle



###########################################################################
#first plot
gradphi=vy
dum     = abs(grad.flatten())
epr_max = max(dum)
	

params = {#'legend.fontsize': 6,
          #'text.latex.preamble': [r"\usepackage{amstext}",],
          #'font',**{'family':'sans-serif','sans-serif':['Helvetica']},
          'axes.linewidth': 2.,
          #'axes.labelsize': 6.5,
          #'text.fontsize': 4,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,	
          'text.usetex': True
          #'figure.figsize': fig_size
          }
mpl.rcParams.update(params)
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
fig1 = plt.figure(1,figsize=(5,5))
ax1 = fig1.add_subplot(111)
epr=ax1.pcolor(grad,cmap='gnuplot2_r',vmin=0,vmax=epr_max)#afmhot_r

epr1=ax1.contour(a,[-0.8,0,0.8],colors='black',alpha = 1.)

ax1.scatter(31.8,68.8,s=3200,fc="None",ec="red",linewidth=2.5,linestyle='solid')

ax1.scatter(52.8,49.8,s=3100,fc="None",ec="green",linewidth=2.5)

ax1.scatter(60.8,20.,s=3100,fc="None",ec="brown",linewidth=2.5)


ax1.scatter(16.7,116.37,s=2800,fc="None",ec="k",linewidth=2.5)

ax1.spines["bottom"].set_linewidth(2)
ax1.spines["left"].set_linewidth(2)
ax1.spines["right"].set_linewidth(2)
ax1.spines["top"].set_linewidth(2)
ax1.set_xticks([])
ax1.set_yticks([])

ax1.plot(26,74.,'o',mfc='c',mec='black',markersize=10,mew=1.5)
ax1.plot(37.3,64.5,'o',mfc='c',mec='black',markersize=10,mew=1.5)

ax1.plot(52.16,45.17,'o',mfc='c',mec='black',markersize=10,mew=1.5)

ax1.plot(51.5,59.6,'^',mfc='r',mec='black',markersize=10,mew=1.5)

ax1.plot(61.89,20.89,'o',mfc='c',mec='black',markersize=10,mew=1.5)

ax1.plot(16.7,118.37,'^',mfc='r',mec='black',markersize=10,mew=1.5)
#ax.annotate('',xy=(30.5,78),xytext=(26,74),color='r',arrowprops=dict(width=2.5,facecolor='r',edgecolor='r', shrink=0.05),zorder=1e3)

#ax.annotate('',xy=(33.5,61.2),xytext=(37.3,64.5),color='r',arrowprops=dict(width=2.5,facecolor='r',edgecolor='r', shrink=0.05),zorder=1e3)

#ax.annotate('',xy=(52.55,40.17),xytext=(52.16,45.0),color='r',arrowprops=dict(width=2.5,facecolor='r',edgecolor='r', shrink=0.05),zorder=1e3)

xp=[52.16,51.5]
yp=[45.17,59.69]


plt.text(-18.50,140.1,r"$(a)$",fontsize=27)

eppr=grad.sum()
ax1.text(-5.5e-1, -9.95, r'$\int \tilde{\sigma} d {\bf r}=%.2g$' %eppr, size=2.5*1e1)
#ax1.text(791.5e-1, -8.95, r'$ \sigma = \psi \nabla ^{4} \psi$', size=2.5*1e1) #( \nabla ^{2} \psi)^{2}
pos = ax1.get_position().get_points()
cax = fig1.add_axes([
     pos[1,1]+0.04, pos[0,1], 0.025,pos[1,1]-pos[0,1] 
])
cbar = fig1.colorbar(epr, cax=cax, orientation='vertical')
cbar.locator = matplotlib.ticker.FixedLocator([0,epr_max])
cbar.ax.yaxis.set_minor_locator(MultipleLocator(epr_max/2))
cbar.update_ticks()
cbar.ax.xaxis.set_ticks_position('top')
#cbar.ax.set_xticklabels(['-1e-7','0','1e-7'])
cbar.ax.set_yticklabels([ r"$0$",r"$0.022$"],fontsize=17.5)
#cbar.ax.set_xticklabels(["-{}".format(fmt(epr_max)),'0',"{}".format(fmt(epr_max))],ha='left')
cbar.ax.set_xlabel(r"${\tilde{\sigma} }(\mathbf{r})$",fontsize=30, rotation=360)
cbar.ax.xaxis.set_label_coords(-0.08,1.13)
#ax1.add_patch(Rectangle((20, 75), 25, 25,edgecolor='black',
#                    facecolor='none',
#                    lw=1.5))



ax1.set_aspect('equal')


plt.savefig('epr_form1.png',bbox_inches='tight',dpi=500)
###########################################################################




fig10 = plt.figure(10,figsize=(5,5))

import matplotlib.patches as mpatches
import matplotlib.path as mpath
grad1=grad#[70:90,14:28]
grad2=grad[55:85,19:47]
ax = fig10.add_subplot(111)
theta = np.linspace(0, 2*np.pi, 400)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = matplotlib.path.Path(verts * radius + center)

#patch = mpatches.PathPatch(path, facecolor='none', edgecolor='k')

vorticity1=ax.pcolor(grad1,cmap='viridis',vmin=0,vmax=0.022,clip_path=(circle, ax.transAxes))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=20., color='white',
        width=0.005, headwidth=0.8, headlength=0.8,clip_path=(circle, ax.transAxes))
#plt.streamplot(X,Y,-vy,-vx, color='white',arrowsize=1,arrowstyle='-|>',density=1)
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["top"].set_linewidth(2)
ax.set_xticks([])
ax.set_yticks([])

from matplotlib.patches import Circle
plt.text(20.51,85.1,r"$(d)$",fontsize=35)

plt.text(25.51,50.1,r"$Simulation$",fontsize=35)

circle=Circle((33,70),14.8,color='r',fill=False,linewidth=3,linestyle='solid')
ax.add_patch(circle)


ax.set_xlabel("ijii")
# pos = ax.get_position().get_points()
# cax = fig10.add_axes([
#     pos[1,0]+0.03, pos[0,0], 0.03,pos[1,1]-pos[0,1]
# ])
# cbar = fig10.colorbar(vorticity1, cax=cax, orientation='vertical')
# cbar.locator = matplotlib.ticker.FixedLocator([0,0.001,0.022])
# cbar.ax.yaxis.set_minor_locator(MultipleLocator(1/2))
# cbar.update_ticks()
# cbar.ax.yaxis.set_ticks_position('right')
# cbar.ax.set_yticklabels([ r"$0$","0.001",r"$0.018$"],fontsize=25)
# #cbar.ax.set_xticklabels([r"$-3 \times 10^{-5}$", r"$0$",r"$3 \times 10^{-5}$"])
# #cbar.ax.set_xticklabels(["-{}".format(fmt(epr_max)),'0',"{}".format(fmt(epr_max))],ha='left')
# cbar.ax.set_xlabel(r"${\sigma}$",fontsize=35, rotation=360) #(\mathbf{r})
# cbar.ax.xaxis.set_label_coords(-0.08,1.13)

ax.annotate('',xy=(30.5,78),xytext=(26,74),color='c',arrowprops=dict(width=2.5,facecolor='c',edgecolor='c', shrink=0.05),zorder=1e3)

ax.annotate('',xy=(33.5,61.2),xytext=(37.3,64.5),color='c',arrowprops=dict(width=2.5,facecolor='c',edgecolor='c', shrink=0.05),zorder=1e3)

#ax.annotate('',xy=(45.5,80.2),xytext=(31.3,69.1),color='r',arrowprops=dict(width=2.5,facecolor='g',edgecolor='g', shrink=0.05),zorder=1e3)


xp=[26,37.3]
yp=[74,64.5]
ax.plot(xp,yp,'--',c='c',linewidth=2.5,zorder=1e3)
ax.scatter(xp,yp,c='r',s=40)

ax.plot(26,74,'o',mfc='c',mec='black',markersize=12,mew=2.)
ax.plot(37.3,64.5,'o',mfc='c',mec='black',markersize=12,mew=2.)

d=np.sqrt((xp[0]-xp[1])**2+(yp[0]-yp[1])**2)
print("d=",d)

xp1=[31.5,43]
yp1=[69.3,80.5]
ax.plot(xp1,yp1,linestyle='dashed',c='r',linewidth=2.5,zorder=1e3)
ax.scatter(xp1,yp1,c='r',s=40)


ax.set_xlim(18,48)#(112,127)#3  
ax.set_ylim(55,85)#(85,100)
ax.axis('off')

R=np.sqrt((xp1[0]-xp1[1])**2+(yp1[0]-yp1[1])**2)
print("R=",R)

print("d/R=",d/R)


import math

def dot(vA, vB):
    return vA[0]*vB[0]+vA[1]*vB[1]

def ang(lineA, lineB):
    # Get nicer vector form
    vA = [(lineA[0][0]-lineA[1][0]), (lineA[0][1]-lineA[1][1])]
    vB = [(lineB[0][0]-lineB[1][0]), (lineB[0][1]-lineB[1][1])]
    # Get dot prod
    dot_prod = dot(vA, vB)
    # Get magnitudes
    magA = dot(vA, vA)**0.5
    magB = dot(vB, vB)**0.5
    # Get cosine value
    cos_ = dot_prod/magA/magB
    # Get angle in radians and then convert to degrees
    angle = math.acos(dot_prod/magB/magA)
    # Basically doing angle <- angle mod 360
    ang_deg = math.degrees(angle)%360
    
    if ang_deg-180>=0:
        # As in if statement
        return 360 - ang_deg
    else: 
        
        return ang_deg


va=30.5-26
vb=78-74

va1=xp[1]-xp[0]
vb1=yp[1]-yp[0]
#linea=[va,vb]
angl=np.arctan(vb/va)
angl1=np.arctan(vb1/va1)
angle=np.arctan(abs((angl1-angl)/(1+angl*angl1)))
print(angl,angl1,angl-angl1)

plt.savefig('ed3.png',bbox_inches='tight',dpi=500)




###########################################################################




fig11 = plt.figure(11,figsize=(5,5))


grad1=grad#[70:90,14:28]
grad2=grad[35:67,35:67]
ax = fig11.add_subplot(111)
theta = np.linspace(0, 2*np.pi, 400)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = matplotlib.path.Path(verts * radius + center)


vorticity1=ax.pcolor(grad1,cmap='viridis',vmin=0,vmax=0.005,clip_path=(circle, ax.transAxes))
plt.quiver(X,Y,np.cos(w),np.sin(w),scale=25., color='white',
        width=0.005, headwidth=0.8, headlength=0.8,clip_path=(circle, ax.transAxes))

        
#plt.streamplot(X,Y,-vy,-vx, color='white',arrowsize=1,arrowstyle='-|>',density=1)
ax.spines["bottom"].set_linewidth(2)
ax.spines["left"].set_linewidth(2)
ax.spines["right"].set_linewidth(2)
ax.spines["top"].set_linewidth(2)
ax.set_xticks([])
ax.set_yticks([])

plt.text(38.51,69.1,r"$(e)$",fontsize=35)

plt.text(43.51,29.5,r"$Simulation$",fontsize=35)

circle=Circle((51.,51.),15.8,color='g',fill=False,linewidth=3,linestyle='solid')
ax.add_patch(circle)

print(grad2.max())

ax.set_xlabel("ijii")
# pos = ax.get_position().get_points()
# cax = fig11.add_axes([
#     pos[1,0]+0.03, pos[0,0], 0.03,pos[1,1]-pos[0,1]
# ])
# cbar = fig11.colorbar(vorticity1, cax=cax, orientation='vertical')
# cbar.locator = matplotlib.ticker.FixedLocator([0,grad2.max()])
# cbar.ax.yaxis.set_minor_locator(MultipleLocator(1/2))
# cbar.update_ticks()
# cbar.ax.yaxis.set_ticks_position('right')
# cbar.ax.set_yticklabels([ r"$0$",r"$0.005$"],fontsize=25)
# #cbar.ax.set_xticklabels([r"$-3 \times 10^{-5}$", r"$0$",r"$3 \times 10^{-5}$"])
# #cbar.ax.set_xticklabels(["-{}".format(fmt(epr_max)),'0',"{}".format(fmt(epr_max))],ha='left')
# cbar.ax.set_xlabel(r"${\sigma}$",fontsize=35, rotation=360) #(\mathbf{r})
# cbar.ax.xaxis.set_label_coords(-0.08,1.13)

ax.annotate('',xy=(52.55,40.17),xytext=(52.16,45.0),color='c',arrowprops=dict(width=2.5,facecolor='c',edgecolor='c', shrink=0.05),zorder=1e3)

ax.annotate('',xy=(33.5,61.2),xytext=(37.3,64.5),color='r',arrowprops=dict(width=2.5,facecolor='r',edgecolor='r', shrink=0.05),zorder=1e3)

#ax.annotate('',xy=(45.5,80.2),xytext=(31.3,69.1),color='r',arrowprops=dict(width=2.5,facecolor='g',edgecolor='g', shrink=0.05),zorder=1e3)

xp=[52.16,51.5]
yp=[45.17,59.69]
ax.plot(xp,yp,'--',c='c',linewidth=2.5,zorder=1e3)
ax.scatter(xp,yp,c='r',s=40)

ax.plot(52.16,45.17,'o',mfc='c',mec='black',markersize=12,mew=2.)

ax.plot(51.5,59.69,'^',mfc='r',mec='black',markersize=12,mew=2.)


xc=[49.18,51.5]
yc=[57.74,59.69]
ax.plot(xc,yc,c='r',linewidth=2.5,zorder=1e3)

xc1=[54.0,51.5]
yc1=[58.01,59.69]
ax.plot(xc1,yc1,c='r',linewidth=2.5,zorder=1e3)


xc2=[51.48,51.5]
yc2=[63.29,59.69]
ax.plot(xc2,yc2,c='r',linewidth=2.5,zorder=1e3)


xp1=[51.81,64.096]
yp1=[51.71,59.69]
ax.plot(xp1,yp1,linestyle='dashed',c='r',linewidth=2.5,zorder=1e3)
ax.scatter(xp1,yp1,c='r',s=40)

ax.set_xlim(35,67)#(112,127)#3  
ax.set_ylim(35,67)#(85,100)
ax.axis('off')

d=np.sqrt((xp[0]-xp[1])**2+(yp[0]-yp[1])**2)
print("d=",d)



R=np.sqrt((xp1[0]-xp1[1])**2+(yp1[0]-yp1[1])**2)
print("R=",R)

print("d/R=",d/R)

va=52.16-52.55
vb=45-40.17

va1=xp[1]-xp[0]
vb1=yp[1]-yp[0]
#linea=[va,vb]
angl=np.arctan(vb/va)
angl1=np.arctan(vb1/va1)
angle=np.arctan(abs((angl1-angl)/(1+angl*angl1)))
print(angl,angl1,angl+angl1)

plt.savefig('ed4.png',bbox_inches='tight',dpi=500)




plt.show()






    
