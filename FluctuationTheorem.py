#!/usr/bin/python

import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import seaborn as sns
import os
from pathlib import Path
#import math
#functions related to elasticity of WLC
def f1(l,lp, kbt):
  return kbt/lp*(1./(4*(1.-l)*(1.-l))-.25+l)

def derf1(l,lp,kbt):
  return kbt/lp*(.5/((1.-l)*(1.-l)*(1.-l))+1.)

def xdefWLC(kbt, l, p, f):
    l0=.9999
    lnew=l0-(f1(l0,p,kbt)-f)/derf1(l0,p,kbt)

    if abs(f)<1.e-5: return 0.0
    while abs(l0-lnew)>1.e-5:
        l0=lnew
        lnew=l0-(f1(l0,p,kbt)-f)/derf1(l0,p,kbt)
    
    xdef=l*lnew        
    return xdef

def intfdexWLC (kbt, l, p, f):
    l0=xdefWLC(kbt,l,p,f)
    return kbt*l/(4.*p)*(1./(1.-l0/l)+2*l0*l0/(l*l)-l0/l)

def stretching_energy(K,l,p,f,fmin):
    total= intfdexWLC(K,l,p,f)-intfdexWLC(K,l,p,fmin)
    return total

forces = [ 5 ,6 ,7 , 8, 9, 10, 11, 12]
delimiter = "   "


p = 0.75   
K = 4.114
l = 0.66
 
#value of mu(pN/nm), N, etc
Deltamu=2.0*K #chemical potential of 1 ATP 3*K
Deltabp=2.0*K #free-energy of 1 bp opening
n_atp=1. #Number of bp opened by 1-step of the motor

force = forces[7]

#frequency bps
k0=1000000.
#montecarlo step = 1e-7
dt=0.5/(k0*5)
fmin=0.0

'''
print (DeltaG)
print(n_atp*Deltamu)
print(Deltabp)
print(2.*force*xdefWLC(K,l,p,force))
print(2.*stretching_energy(K,l,p,force,fmin))
#time.sleep(30.5)
print("pre")
print(k0*np.exp(-DeltaG/(2.*K))*dt)
print("pun")
print(k0*np.exp(DeltaG/(2.*K))*dt)
'''



for i in range(1,2):
	force=forces[i]
	fmin=0.0
	DeltaG=n_atp*Deltamu-Deltabp+2.*force*xdefWLC(K,l,p,force)-2.*stretching_energy(K,l,p,force,fmin)
	Ateo=(n_atp*Deltamu/(2*xdefWLC(K,l,p,force))-Deltabp/(2*xdefWLC(K,l,p,force))+force-1.*stretching_energy(K,l,p,force,fmin)/xdefWLC(K,l,p,force))/K
	print(Ateo)
	for j in range (0,1):
		#CHECK parameter dt is the same than for the files processed in hist
		Dt=10*dt*2**j
		#print(Dt)
		dn=10*2**j
		n=0
		folder='/home/xavi/Documents/Helicase/Simulations/Simulation/Simulation_221018_4/FT/DATA/'
		filein='%s%s_%dhist.txt' %(folder,force,dn)
		#File_in Dx, p(Dx)
		foldout='/home/xavi/Documents/Helicase/Simulations/Simulation/Simulation_221018_4/FT/'
		#reading input file
		data=np.loadtxt(filein)
		data2=np.asarray(data)
		dx=data2[:,0]
		pdx=data2[:,1]
		
		#computing the average <Dx>
		dxav=0.
		pdxav=0.
		##Cal renormalitzar de nou les probabilitats (al descartar el 0 la integral ha canviat)
		for k in range(0,len(dx)):
			dxav+=dx[k]*pdx[k]
			pdxav+=pdx[k]
		dxav/=pdxav
		pdx/=pdxav
		
		print(force)
		print(Dt)
		print(dxav)
		print("\n")
		
		epsilon=0.08
		delta=0.4
		y=[]
		x=[]
		for k in range(0,len(dx)):
			if (dx[k]<=0.) & (dx[k]<=-delta):
				for l in range(0,len(dx)):
				#comprovem si x[l]+x[k]<epsilon --> computem 
					if (dx[l]>=0.) & (dx[l]>=delta):
						if np.absolute(dx[k]+dx[l]) <= epsilon:
							x.append(dx[l]/dxav)
							y.append((np.log(pdx[l]/pdx[k]))/dxav)
		'''		
		print(x)
		#print("\n")
		#print(y)
		print(len(x))
		print("\n")
		print(len(y))
		#plt.plot(x,y, 'o')
		plt.plot(x,y)
		plt.show()
		'''
		#writing to an output file for each DT		
		fileout='%s%s_%dFT.txt' %(foldout,force,dn)
		fout=open(fileout,"w")
		#fout.write()
		for k in range(0,len(x)):
			fout.write('%.12f %.12f\n'%(x[k],y[k]))		
		fout.close()
