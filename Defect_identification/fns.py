###############################################################
import numpy as np




###############################################################
def gaussian(distance,dx,cutoff):
	#x=lattice_angle.shape[0]
	#y=lattice_angle.shape[1]
	
	weight =np.exp(-(distance**2)/(2*cutoff**2*dx**2))
	normalization=np.sqrt(2*np.pi*cutoff**2*dx**2)
	gaussian=weight/normalization
	
	
	return gaussian
Nx=128	
def next(i,j):
	return (i+j) % Nx
	
def back(i,j):

	if (i>j-1):
		return (i-j)
	else:
		return Nx+(i-j)
	
def smoothening(s11,s12,Lx,dx,cutoff):
	
	#calculating positions
	x=np.arange(0,Lx,dx)
	y=np.arange(0,Lx,dx)
	#X,Y=np.meshgrid(x,y)
	Nx=int(Lx/dx)
	#X.reshape(Nx*Nx)
	#Y.reshape(Nx*Nx)

	xlist=[]
	ylist=[]
	for j in range(Nx):
		for i in range(Nx):
			xlist.append(x[i])
			ylist.append(y[j])
		
	sm_s11=np.zeros((Nx,Nx))
	sm_s12=np.zeros((Nx,Nx))	
			
	#
	for particle in range(0,int((Lx/dx)**2)):
		#for i in range()
		jx=int((xlist[particle]/dx) % Nx)
		jy=int((ylist[particle]/dx) % Nx)
		
		for djy in range(0,cutoff):
		
			djx=0
			while(djx*djx+djy*djy < cutoff ):
			
				distance=(djx*djx+djy*djy)*dx**2
				
				#firts quadrant
				sm_s11[next(jx,djx),next(jy,djy)]   +=s11[jx,jy] * gaussian(distance,dx,cutoff)
				sm_s12[next(jx,djx),next(jy,djy)]   +=s12[jx,jy] * gaussian(distance,dx,cutoff)
				
				#third quadrant
				sm_s11[back(jx,djx),back(jy,djy)]   +=s11[jx,jy] * gaussian(distance,dx,cutoff)
				sm_s12[back(jx,djx),back(jy,djy)]   +=s12[jx,jy] * gaussian(distance,dx,cutoff)
				
				#fourth quadrant
				sm_s11[next(jx,djx),back(jy,djy)]   +=s11[jx,jy] * gaussian(distance,dx,cutoff)
				sm_s12[next(jx,djx),back(jy,djy)]   +=s12[jx,jy] * gaussian(distance,dx,cutoff)
				
				#second quadranr
				sm_s11[back(jx,djx),next(jy,djy)]   +=s11[jx,jy] * gaussian(distance,dx,cutoff)
				sm_s12[back(jx,djx),next(jy,djy)]   +=s12[jx,jy] * gaussian(distance,dx,cutoff)
				
				djx +=1
				
	lattice_angle =np.arctan2(sm_s11,sm_s12)*.5#+np.pi/4
				
							
	return lattice_angle
	
def smoothening_vort(v,Lx,dx,cutoff):
	
	#calculating positions
	x=np.arange(0,Lx,dx)
	y=np.arange(0,Lx,dx)
	#X,Y=np.meshgrid(x,y)
	Nx=int(Lx/dx)
	#X.reshape(Nx*Nx)
	#Y.reshape(Nx*Nx)

	xlist=[]
	ylist=[]
	for j in range(Nx):
		for i in range(Nx):
			xlist.append(x[i])
			ylist.append(y[j])
		
	sm_v=np.zeros((Nx,Nx))
	#sm_s12=np.zeros((Nx,Nx))	
			
	#
	for particle in range(0,int((Lx/dx)**2)):
		#for i in range()
		jx=int((xlist[particle]/dx) % Nx)
		jy=int((ylist[particle]/dx) % Nx)
		
		for djy in range(0,cutoff):
		
			djx=0
			while(djx*djx+djy*djy < cutoff ):
			
				distance=(djx*djx+djy*djy)*dx**2
				
				#firts quadrant
				sm_v[next(jx,djx),next(jy,djy)]   +=v[jx,jy] * gaussian(distance,dx,cutoff)
				#sm_s12[next(jx,djx),next(jy,djy)]   +=s12[jx,jy] * gaussian(distance,dx,cutoff)
				
				#third quadrant
				sm_v[back(jx,djx),back(jy,djy)]   +=v[jx,jy] * gaussian(distance,dx,cutoff)
				#sm_s12[back(jx,djx),back(jy,djy)]   +=s12[jx,jy] * gaussian(distance,dx,cutoff)
				
				#fourth quadrant
				sm_v[next(jx,djx),back(jy,djy)]   +=v[jx,jy] * gaussian(distance,dx,cutoff)
				#sm_s12[next(jx,djx),back(jy,djy)]   +=s12[jx,jy] * gaussian(distance,dx,cutoff)
				
				#second quadranr
				sm_v[back(jx,djx),next(jy,djy)]   +=v[jx,jy] * gaussian(distance,dx,cutoff)
				#sm_s12[back(jx,djx),next(jy,djy)]   +=s12[jx,jy] * gaussian(distance,dx,cutoff)
				
				djx +=1
				
	#lattice_angle =np.arctan2(sm_s11,sm_s12)*.5+np.pi/4
				
							
	return sm_v
				
	




def mod_angle(x,y):

	if np.abs(x-y)< np.pi/2:
		return x-y
	else:
		return (x-y-np.sign(x-y)*np.pi)#*2*np.pi



def lattice_charge(lattice_angle,threshold):
	
	
	lattice_charge=np.zeros((Nx,Nx))
	dx=0.5
	xvec=[]
	yvec=[]
	cvec=[]
	for i in range(0,Nx):
		x=(dx*(2*i+1)/2)
		xvec.append(x)
		for j in range(0,Nx):
			y=(dx*(2*j+1)/2)
			yvec.append(y)
			
			charge  =mod_angle(lattice_angle[back(i,1),j],lattice_angle[back(i,1),next(j,1)])  #left -> left-top
			charge += mod_angle(lattice_angle[back(i,1),next(j,1)],lattice_angle[i,next(j,1)]) #left-top -> top
			charge += mod_angle(lattice_angle[i,next(j,1)],lattice_angle[next(i,1),next(j,1)]) #top -> top-right
			charge += mod_angle(lattice_angle[next(i,1),next(j,1)],lattice_angle[next(i,1),j]) #top-right -> right
			charge += mod_angle(lattice_angle[next(i,1),j],lattice_angle[next(i,1),back(j,1)]) #right -> right-bottom
			charge += mod_angle(lattice_angle[next(i,1),back(j,1)],lattice_angle[i,back(j,1)]) #right-bottom -> bottom
			charge += mod_angle(lattice_angle[i,back(j,1)],lattice_angle[back(i,1),back(j,1)]) #bottom -> bottom-left
			charge += mod_angle(lattice_angle[back(i,1),back(j,1)],lattice_angle[back(i,1),j]) #bottom-left -> left
			
			
			if abs(charge) > threshold:
				
				lattice_charge[i,j]=charge/(2*np.pi)
				cvec.append(charge)
				
	return lattice_charge #xvec,yvec,cvec #
			
			



















	
