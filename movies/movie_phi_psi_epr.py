from pylab import *
from matplotlib.ticker import  *
from matplotlib.pyplot import  *
from matplotlib import animation
from matplotlib import cm
from scipy.io import *
import time
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as mticker




# References
# ----------
# https://matplotlib.org/examples/animation/rain.html
# https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/

fig_width_pt  = 200. #246
inches_per_pt = 1./72
fig_width 	  = fig_width_pt*inches_per_pt
fig_size	 	  = [3*fig_width, fig_width]

params = {'legend.fontsize': 6,
          #'text.latex.preamble': [r"\usepackage{amstext}",],
          'axes.linewidth': 1.5,
          'axes.labelsize': 6.5,
          #'text.fontsize': 4,
          'xtick.labelsize': 6,
          'ytick.labelsize': 6,	
          'text.usetex': True,
          'figure.figsize': fig_size}
rcParams.update(params)


left   = 3e-2
bottom = 1.4e-1
width  = 2.8e-1
height = 8.4e-1

leftb   = 3.6e-1
bottomb = bottom
widthb  = width
heightb = height


leftc   = 6.9e-1
bottomc = bottom
widthc  = width
heightc = height

f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.1e' % x))
fmt = mticker.FuncFormatter(g)

dx=0.5000
# compute the energy spectrum numerically
def energy_spectrum(nx,ny,psi2,phi2):
    
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
    
    kx[0:int(nx/2)] = np.float64(np.arange(0,int(nx/2))) *2*np.pi/(np.float64(nx)*dx)
    kx[int(nx/2):nx] = np.float64(np.arange(-int(nx/2),0))* 2*np.pi/(np.float64(nx)*dx)

    ky[0:ny] = kx[0:ny]
    
    kx[0] = epsilon
    ky[0] = epsilon
    kx, ky = np.meshgrid(kx, ky, indexing='ij')
    
    #a = pyfftw.empty_aligned((nx,ny),dtype= 'complex128')
    #b = pyfftw.empty_aligned((nx,ny),dtype= 'complex128')

    #fft_object = pyfftw.FFTW(a, b, axes = (0,1), direction = 'FFTW_FORWARD')
    
    #wf = fft_object(w[0:nx,0:ny]) 
    #vxf = fft_object(vx[0:nx,0:ny] )
    #vyf = fft_object(vy[0:nx,0:ny]) 
    w=np.zeros((psi2.shape[0], nx*ny)) 
    vxs=np.zeros((psi2.shape[0], nx*ny)) 
    vys=np.zeros((psi2.shape[0], nx*ny)) 
    epr=np.zeros((psi2.shape[0], nx*ny)) 
    time_s=int(psi2.shape[0])
    
    for i in range(psi2.shape[0]):
    	psi=psi2[i,:].reshape(nx,nx)
    	phi=phi2[i,:].reshape(nx,nx)
    	psif=np.fft.fft2(psi)
    	phif=np.fft.fft2(phi)
    
    #print(len(vxf))
    
    
    	kk = np.sqrt(kx[:,:]**2 + ky[:,:]**2)
    	kx2=kx[:,:]**2
    	ky2=ky[:,:]**2
    #s11=0.5*(- ky2*psif + kx2 * psif)
    #s12=  -(1j*kx[:,:] * psif) * (1j*ky[:,:]* psif)
    	d2psi_ft= - kk**2 *psif
    	d2phi_ft= - kk**4 *phif

    	vx_ft= 1j*ky[:,:]* psif
    	vy_ft= -1j*kx[:,:]* psif
    	
    	dypsi_ft=(1j*ky[:,:]* psif)
    	dxpsi_ft=(1j*kx[:,:]* psif)

    	d2psi=np.real(np.fft.ifft2(d2psi_ft))#/(0.003/1.67)
    	d2phi=np.real(np.fft.ifft2(d2phi_ft))#/(0.003/1.67)

    	vx=np.real(np.fft.ifft2(vx_ft))
    	vy=np.real(np.fft.ifft2(vy_ft))
    	
    	w[i,:]=d2psi.reshape(nx*ny)
    	vxs[i,:]=vx.reshape(nx*ny)
    	vys[i,:]=vy.reshape(nx*ny)
    	
    	dxpsi=np.real(np.fft.ifft2(dxpsi_ft))

    	dypsi=np.real(np.fft.ifft2(dypsi_ft))


    	#S11=0.5 * ((dypsi **2)-(dxpsi)**2)
    	#S12= - dxpsi * dypsi
    
    	#r=S11**2+S12**2
    
    	#s11_ft=np.fft.fft2(S11)
    	#s12_ft=np.fft.fft2(S12)
    
    	#dxs11_ft=(1j*kx[:,:]* s11_ft)
    	#dys11_ft=(1j*ky[:,:]* s11_ft)
    
    	#dxs12_ft=(1j*kx[:,:]* s12_ft)
    	#dys12_ft=(1j*ky[:,:]* s12_ft)
    
    	#dxs11=np.real(np.fft.ifft2(dxs11_ft))
    	#dys11=np.real(np.fft.ifft2(dys11_ft))
    	#dxs12=np.real(np.fft.ifft2(dxs12_ft))
    	#dys12=np.real(np.fft.ifft2(dys12_ft))
    
    	#rho=dxs11*dys12-dxs12-dys11
    
    
    
    #grad term
    
    	d2dxpsi_ft= - kk **2 *dxpsi_ft
    	d2dypsi_ft= - kk **2 *dypsi_ft
    
    	d2dxpsi=np.fft.ifft2(d2dxpsi_ft)
    	d2dypsi=np.fft.ifft2(d2dypsi_ft)
    
    	grad=d2phi**2#*1.67*.103/(.003) 
    	epr[i,:]=grad.reshape(nx*ny)
    


    
    
        
    return w,epr,time_s



# =======================
#  Load data
# =======================

idx       = [0]
#param_all = ['N_64_T_1e-5_con']
#param_all = ['N_64_T_5e-3_ext']
param_all = ['lam']#Turb_k1_.05_movie

N_plot = len(param_all)

# =======================
#  Plot data
# =======================

# Start time counter
start = time.time()

for j in range(N_plot):
	param  = param_all[j]

	# Load parameters
	d_name          = param	
	data            = loadmat(d_name + '.mat')
	phi             = data['data_phi']
	psi             = data['data_psi']
	#psi             = data['data_ome']
	#epr             = data['data_psi']
	N               = int(data['Nx'])
	k1             = data['kappa1']
	a               = data['a']
	b               = data['b']
	k               = data['kappa']
	eta               = data['eta']
	#totaltime       = data['totaltime']
	#interval_record = data['interval_record']
	#nb_record       = int(totaltime/interval_record)

	# Plot parameters
	#if a<0:
	#	phi_max = sqrt(-a/b)
	
	print("============================== \n Activity (k')=%0.3f \n a =%0.3f \n b =%0.3f \n k =%0.3f \n eta =%0.3f" %(k1,a,b,k,eta))
	print("============================== \n Making movie")
	phi_max = 1
	psi_norm=abs(float(k1))/float(eta)
	w,grad,time_s=energy_spectrum(N,N,phi,psi)
	#grad=grad
	dum     = abs(phi.flatten())
	phi_max = max(dum)
	dum     = abs(psi.flatten())
	psi_max = max(dum)/psi_norm
	dum     = abs(grad.flatten())
	epr_max = max(dum)




	# Plot data: phi
	fig = figure()
	ax  = axes([left, bottom, width, height])
	#cax = axes([left_cb, bottom_cb, width_cb, height_cb])
	ax.set_xticks([])
	ax.set_yticks([])

	density = ax.pcolor(phi[0,:].reshape(N,N), cmap='RdBu', vmin=-phi_max, vmax=phi_max)#RdBu
	#cb      = colorbar(density, cax, orientation='horizontal')
	#cb.ax.tick_params(labelsize=1e0, length=0, pad=1e1)
	props=dict(boxstyle='square',facecolor='white',alpha=1.)
	#ax.text(0.05,0.95,"",transform=ax.transAxes,fontsize=12,verticalalignment='top',bbox=props)
	txt_time = text(0.03,0.03, '', transform=ax.transAxes,fontsize=12,verticalalignment='bottom',bbox=props) #8e-1, 2.5
	
	txt_zeta = text(0.03,0.91, r"$\zeta = %0.3f$"%(float(k1)-float(k)), transform=ax.transAxes,fontsize=12,verticalalignment='bottom',bbox=props) #8e-1, 2.5
	
	#text(-3e-1, 0, r'$\phi=-%.2g$' %phi_max, size=1e1)
	#text(1.05e0, 0, r'$%.2g$' %phi_max, size=1e1)
	pos = ax.get_position().get_points()
	cax = fig.add_axes([
    	pos[0,0], pos[0,0]+0.07, pos[1,0]-pos[0,0], 0.025])
	cbar = colorbar(density, cax=cax, orientation='horizontal')
	cbar.locator = matplotlib.ticker.FixedLocator([-phi_max,0,phi_max])
	cbar.ax.xaxis.set_minor_locator(MultipleLocator(phi_max/2))
	cbar.update_ticks()
	cbar.ax.xaxis.set_ticks_position('bottom')
	#cbar.ax.set_yticklabels(['-1e-7','0','1e-7'])
	cbar.ax.set_xticklabels([r"$-%0.2g$" %phi_max,r"$0$",r"$%0.2g$" %phi_max],fontsize=12)
	cbar.ax.set_ylabel(r"${\phi}$",fontsize=17, rotation=360)
	cbar.ax.yaxis.set_label_coords(-0.05,.1)






	# Plot data: psi
	axb  = axes([leftb, bottomb, widthb, heightb])
	#caxb = axes([left_cbb, bottom_cbb, width_cbb, height_cbb])
	axb.set_xticks([])
	axb.set_yticks([])

	vorticity = axb.pcolor(psi[0,:].reshape(N,N)/psi_norm, cmap='BrBG',vmin=-psi_max, vmax=psi_max)#vmin=-psi_max, vmax=psi_max#BrBG
	#cbb       = colorbar(vorticity, caxb, orientation='horizontal')
	#cbb.ax.tick_params(labelsize=1e0, length=0, pad=1e1)
	#text(-3.5e-1, 0, r'$\psi=-%.2g$' %psi_max, size=1e1)
	#text(1.05e0, 0, r'$%.2g$' %psi_max, size=1e1)
	
	pos = axb.get_position().get_points()
	cbx = fig.add_axes([
    	pos[0,0], pos[0,0]-0.26, pos[1,0]-pos[0,0], 0.025])
	cbar = colorbar(vorticity, cax=cbx, orientation='horizontal')
	cbar.locator = matplotlib.ticker.FixedLocator([-psi_max/1.0,0,psi_max/1.0])
	cbar.ax.xaxis.set_minor_locator(MultipleLocator(psi_max/2))
	cbar.update_ticks()
	cbar.ax.xaxis.set_ticks_position('bottom')
	#cbar.ax.set_yticklabels(["val = -{}".format(fmt(psi_max)),'0',"val = {}".format(fmt(psi_max))])
	cbar.ax.set_xticklabels([r"$- %0.3g$" %psi_max,r"$0$",r"$%0.3g$" %psi_max],ha='center',fontsize=12)#([r"$-3 \times 10^{-5}$", r"$0$",r"$3 \times 10^{-5}$"])
	#cbar.ax.set_xticklabels(["-{}".format(fmt(psi_max/1.0)),'0'," {}".format(fmt(psi_max/1.0))],ha='center')#([r"$-3 \times 10^{-5}$", r"$0$",r"$3 \times 10^{-5}$"])
	cbar.ax.set_ylabel(r"${\tilde{\psi}}$",fontsize=17, rotation=360)
	cbar.ax.yaxis.set_label_coords(-0.05,.1)

	# Plot data: vorticity
	axc  = axes([leftc, bottomc, widthc, heightc])
	#caxc = axes([left_cbc, bottom_cbc, width_cbc, height_cbc])
	axc.set_xticks([])
	axc.set_yticks([])

	vort = axc.pcolor(grad[0,:].reshape(N,N), cmap='gnuplot2_r',vmin=0, vmax=epr_max)#, vmin=-epr_max, vmax=epr_max
	#cbc  = colorbar(vort, caxc, orientation='horizontal')
	#cbc.ax.tick_params(labelsize=1e0, length=0, pad=1e1)
	#text(-3.5e-1, 0, r'$\sigma=-%.2g$' %epr_max, size=1e1)
	#text(1.05e0, 0, r'$%.2g$' %epr_max, size=1e1)
	txt_epr = text(-5e-1, 3e0, '', size=1.5e1)
	
	pos = axc.get_position().get_points()
	ccx = fig.add_axes([
    	pos[0,0], pos[0,0]-0.59, pos[1,0]-pos[0,0], 0.025])
	cbar = colorbar(vort, cax=ccx, orientation='horizontal')
	cbar.locator = matplotlib.ticker.FixedLocator([0,epr_max])
	cbar.ax.xaxis.set_minor_locator(MultipleLocator(epr_max/2))
	cbar.update_ticks()
	cbar.ax.xaxis.set_ticks_position('bottom')
	#cbar.ax.set_yticklabels(["val = -{}".format(fmt(epr_max)),'0',"val = -{}".format(fmt(epr_max))])
	#cbar.ax.set_xticklabels(["-{}".format(fmt(epr_max)),'0',"{}".format(fmt(epr_max))],ha='center')
	cbar.ax.set_xticklabels([r"0",r"$%0.3g$" %epr_max],ha='right',fontsize=12)
	#([r"$-3 \times 10^{-5}$", r"$0$",r"$3 \times 10^{-5}$"])
	#cbar.ax.set_xticklabels(['','',"{}".format(fmt(epr_max))],ha='right')
	cbar.ax.set_ylabel(r"${\dot{\sigma}}$",fontsize=17, rotation=360)
	cbar.ax.yaxis.set_label_coords(-0.05,.1)
	
	
	# Animate plot
	def animate(i):
			txt_time.set_text(r'$t = %.4g \times 10^{2}$' %(3*i))
			txt_epr.set_text(r'$\int\dot{\sigma} d{\bf r} = %.3g$' %sum(grad[i,:]))
			density.set_array(phi[i,:])
			#density.set_clim(-phi[i,:].max(),phi[i,:].max())
			
			
			
		
			vorticity.set_array(psi[i,:]/psi_norm)
			#vorticity.set_clim(-psi[i,:].max(),psi[i,:].max())
			
			vort.set_array(grad[i,:])
			#vort.set_clim(-grad[i,:].max(),grad[i,:].max())
			#heat.set_array(-w[i,:])
			print (d_name + ': %d / %d'%(i,time_s))
			return density, vorticity,vort, txt_epr,txt_time

	ani = animation.FuncAnimation(fig, animate, 600)#, interval=40, blit=True)

	# Save movie	
	movie_name = 'movie_check_' + param + '.mp4'
	ani.save(movie_name, fps=30, dpi=600)#, bitrate=1e3)

# End script message
duration = time.time() - start
print ('Plot duration:', str('%.3g' %duration), 'seconds')


