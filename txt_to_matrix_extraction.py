import sys, os, re, math
from scipy.io import savemat
import numpy as np


# References
# https://docs.python.org/3/library/os.html
# https://docs.python.org/3/library/sys.html
# https://scipy.org/



nstart=0	# Start timestep
nint=500000	# Increment
nend=11000000	# final timestep


#flags to extract data
phi_files=1
psi_files=1
epr_files=1


ntotal=int((nend-nstart)/nint)+1

#reading log.txt file and absorbing parameters
with open("log.txt", "r") as file:
	Nx_b=False
	Ny_b=False
	Lx_b=False
	Ly_b=False
	dt_b=False
	#dx_b=False
	
	#D_b=False
	kappa_b=False
	kappa1_b=False
	lambda_b=False
	m_b=False
	zeta_b=False
	
	eta_b=False
	a_b=False
	b_b=False
	tem_b=False
	
	
	for line in file:
        	if ":" in line:
            		a, b = map(str.strip, line.split(":"))
            		if a=='Nx':
            			Nx=int(b)
            			Nx_b=True
            		elif a=='Ny':
            			Ny=int(b)
            			Ny_b=True
            		elif a=='Lx':
            			Lx=b
            			Lx_b=True
            		elif a=='Ly':
            			Ly=b
            			Ly_b=True
            		elif a=='dt':
            			dt=b
            			dt_b=True
            		elif a=='lambda':
            			Lambda=b
            			lambda_b=True
            		elif a=='kappa':
            			kappa=b
            			kappa_b=True
            		elif a=='kappa1':
            			kappa1=b
            			kappa1_b=True
            		elif a=='m':
            			m=b
            			m_b=True
            		elif a=='zeta':
            			zeta=b
            			zeta_b=True
            		elif a=='eta':
            			eta=b
            			eta_b=True
            		elif a=='a':
            			a1=b
            			a_b=True
            		elif a=='b':
            			b1=b
            			b_b=True
            		elif a=='tem':
            			tem=b
            			tem_b=True
            		
            			if (Nx_b & Ny_b & Lx_b & Ly_b & dt_b & lambda_b & kappa_b & kappa1_b & zeta_b & eta_b & m_b & a_b & b_b & tem_b & eta_b):
            				break


# Set arrays to save data
phi=np.zeros((ntotal, Nx*Ny)) 
psi=np.zeros((ntotal, Nx*Ny)) 
epr=np.zeros((ntotal, Nx*Ny)) 


ii=0
if (phi_files==1):
	for i in range(nstart,nend+nint,nint):
		name='phi_'+str(i)+'.txt'     #name each data file
		if os.path.exists(name):
			c=np.loadtxt(name)
			
			phi[ii,:]=c.reshape(Nx*Ny)
			ii=ii+1

			#with open(s ,'r') as f:
			#	list3=[]
			#	for line in f:
			#		strip_lines=line.strip()
			#		listli=strip_lines.split()
			#		m=list3.append(listli)
			#	d=np.asfarray(list3,float)
		
		else:
			print('file', name, 'does not exists')
			
ii=0
if (psi_files==1):
	for i in range(nstart,nend+nint,nint):
		name='psi_'+str(i)+'.txt'
		if os.path.exists(name):
			c=np.loadtxt(name)
			
			psi[ii,:]=c.reshape(Nx*Ny)
			ii=ii+1

			#with open(s ,'r') as f:
			#	list3=[]
			#	for line in f:
			#		strip_lines=line.strip()
			#		listli=strip_lines.split()
			#		m=list3.append(listli)
			#	d=np.asfarray(list3,float)
		
		else:
			print('file', name, 'does not exists')
		
		
#saving data as a matrix
savemat('Data_matrix_Active_H.mat', {'Nx':Nx,'Ny':Ny,'Ny':Lx,'Ny':Ly,'dt':dt,'lambda':Lambda,'kappa':kappa,'kappa1':kappa1,'zeta':zeta,'eta':eta,'m':m,'a':a1,'b':b1,'temp':tem,'data_phi':phi, 'data_psi':psi, 'data_epr':epr})


print('Completed succesfully')
