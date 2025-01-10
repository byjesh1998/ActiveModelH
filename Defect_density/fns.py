###############################################################
import numpy as np
from scipy.spatial import distance




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
				
	lattice_angle =np.arctan2(sm_s11,sm_s12)*.5
				
							
	return lattice_angle
				
	




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
	pd=[]
	nd=[]
	for i in range(0,Nx):
		#x=(dx*(2*i+1)/2)
		#xvec.append(x)
		for j in range(0,Nx):
			
			
			charge  =mod_angle(lattice_angle[back(i,1),j],lattice_angle[back(i,1),next(j,1)])  #left -> left-top
			charge += mod_angle(lattice_angle[back(i,1),next(j,1)],lattice_angle[i,next(j,1)]) #left-top -> top
			charge += mod_angle(lattice_angle[i,next(j,1)],lattice_angle[next(i,1),next(j,1)]) #top -> top-right
			charge += mod_angle(lattice_angle[next(i,1),next(j,1)],lattice_angle[next(i,1),j]) #top-right -> right
			charge += mod_angle(lattice_angle[next(i,1),j],lattice_angle[next(i,1),back(j,1)]) #right -> right-bottom
			charge += mod_angle(lattice_angle[next(i,1),back(j,1)],lattice_angle[i,back(j,1)]) #right-bottom -> bottom
			charge += mod_angle(lattice_angle[i,back(j,1)],lattice_angle[back(i,1),back(j,1)]) #bottom -> bottom-left
			charge += mod_angle(lattice_angle[back(i,1),back(j,1)],lattice_angle[back(i,1),j]) #bottom-left -> left
			
			
			if abs(charge) > threshold:
				
				charge2  =mod_angle(lattice_angle[i,j],lattice_angle[i,next(j,1)])
				charge2 +=mod_angle(lattice_angle[i,next(j,1)],lattice_angle[next(i,1),next(j,1)])
				charge2 +=mod_angle(lattice_angle[next(i,1),next(j,1)],lattice_angle[next(i,1),j])
				charge2 +=mod_angle(lattice_angle[next(i,1),j],lattice_angle[i,j])
				
				if abs(charge2) > threshold:
					if charge2 < 0:
					
						lattice_charge[i,j]=charge/(2*np.pi)
						x=((2*j+1)/2)#(dx*(2*i+1)/2)
						y=((2*i+1)/2)#(dx*(2*j+1)/2)
						pd.append((x,y))
					else:
						#lattice_charge[i,j]=charge/(2*np.pi)
						x=((2*j+1)/2)#(dx*(2*i+1)/2)
						y=((2*i+1)/2)#(dx*(2*j+1)/2)
						nd.append((x,y))
					
					#yvec.append(i)
			
					#x=(dx*(2*i+1)/2)
					#xvec.append(j)
					#cvec.append(charge/(2*np.pi))
				
	return pd,nd ,lattice_charge
			
			


def defect_anhilation(pd,nd,minimum_dist):


	defect_dist=distance.cdist(pd,nd)
	min_dist=defect_dist.min()
	min_dist_ind=np.where(defect_dist == np.min(defect_dist))
	
	pd=np.delete(pd,min_dist_ind[0],axis=0)
	nd=np.delete(nd,min_dist_ind[1],axis=0)
	
	while(min_dist < minimum_dist):
	
		defect_dist=distance.cdist(pd,nd)
		min_dist=defect_dist.min()
		min_dist_ind=np.where(defect_dist == np.min(defect_dist))
	
		pd=np.delete(pd,min_dist_ind[0],axis=0)
		nd=np.delete(nd,min_dist_ind[1],axis=0)
		
	return pd,nd
















	
