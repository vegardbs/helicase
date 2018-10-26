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
Deltamu=3.0*K #chemical potential of 1 ATP 3*K
Deltabp=2.0*K #free-energy of 1 bp opening
n_atp=1. #Number of bp opened by 1-step of the motor

force = forces[7]

#frequency bps
k0=1000000.
#montecarlo step = 1e-7
dt=0.5/(k0*5)
fmin=0.0

DeltaG=n_atp*Deltamu-Deltabp+2.*force*xdefWLC(K,l,p,force)-2.*stretching_energy(K,l,p,force,fmin)
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



for i in range(0,8):
	force=forces[i]
	folder='../%s/*.dat'%(force)
	#average file
	averagefile='../%s/average.txt'%(force)
	var_averagefile='../%s/var_average.txt'%(force)
	my_file=Path(averagefile)
	cont=0
	numfiles=0
	foldout='/home/xavi/Documents/Helicase/Simulations/Simulation/Simulation_221018_4/'
	for filename in glob.glob(folder):
			g=open(filename,"r")
			lines=g.readlines()
			for y in lines:
				cont+=1
			g.close()
			numfiles+=1
			if numfiles>=1:
				break;
	numlines=cont
	numfilestotal=numfiles
	print("Numfile: ",numfilestotal)
	print("Numlines: ", numlines)
	histDt=[8]
	for j in range (0,6):
		Dt=10*dt*2**j
		#print(Dt)
		dn=10*2**j
		n=0
		#print(dn)
		#for k in range(1,numfilestotal)
		for k in range(0,1):
			print(forces[i])
			print(k)
			data=np.loadtxt('/home/xavi/Documents/Helicase/Simulations/Simulation/Simulation_221018_4/' + str(forces[i]) + '/'+ str(k)+'.dat')
			data2 = np.asarray(data)
			x = data2[:,1]
			t = data2[:,0]
			n1=0
			n2=dn
			
			N_steps =int(len(t)-dn)
			dbp = np.zeros(N_steps) 
			while n2 < N_steps:
				dbp[n] = x[int(n2)] - x[int(n1)]
				n1 += 1
				n2 += 1
				n += 1

		outputfile='%s%s_%d.txt'%(foldout,force,dn)
		fout=open(outputfile,"w")
		for k in range(0,n):		
			fout.write('%.12f \n'%(dbp[k]))
		fout.close()
		#dbp_clean=np
		nbin=500
		xmin=-20.0
		xmax=40.0
		hist, bin_edges = np.histogram(dbp,bins=nbin, range=(xmin,xmax), density=True)
		bin_width = bin_edges[1]-bin_edges[0]
		bin_center = bin_edges + bin_width/2
		bin_center = bin_center[0:-1]
		x=np.zeros(nbin)
		#print (bin_edges[0])
		#print(bin_edges)
		for n in range (0,nbin):
    			x[n]=bin_edges[0]+bin_width*n

		#plt.plot(x,hist)
		#plt.show()
		cleanhist=hist[hist!=0]
		cleanx=[]
		for n in range (0,nbin):
			if hist[n]!=0:
				cleanx.append(x[n])
		#cleanx=x[]
		#print(hist)
		#print(cleanhist)
		#print(x)
		#print(cleanx)
		#plt.plot(cleanx,cleanhist)
		#plt.show()
		outputfilehist='%s%s_%dhist.txt'%(foldout,force,dn)
		fout=open(outputfilehist,"w")
		for k in range(0,len(cleanx)):		
			fout.write('%.12f %.12f\n'%(cleanx[k],cleanhist[k]))
		fout.close()
