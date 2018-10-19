
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 12:14:52 2018
@author: vegard
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

# Experimental parameters

p = 0.75   
K = 4.114
l = 0.66
forces = [ 5. ,6. ,7. , 8., 9., 10., 11., 12.]
gammas = [ 0.5538,    0.6155,    0.6668  ,  0.7100 ,   0.7466   , 0.7781  ,  0.8054 ,  0.8293]   
gammas_phi = [  0.3532,    0.4724,    0.5741,    0.6537,    0.7139,    0.7596,    0.7951,    0.8293]

 
#value of mu(pN/nm), N, etc
Deltamu=20.0*K #chemical potential of 1 ATP
Deltabp=2.0*K #free-energy of 1 bp opening
n_atp=1. #Number of bp opened by 1-step of the motor

force = forces[6]
gamma = gammas[0]

#frequency bps
k0=1000000.
#montecarlo step
dt=0.5/k0
fmin=0.0
   
N_steps = 100000
N_traces = 1

#print(N*Deltamu-Deltabp+2*force*xdefWLC(K,l,p,force)-2*stretching_energy(K,l,p,force,fmin))
DeltaG=n_atp*Deltamu-Deltabp+2.*force*xdefWLC(K,l,p,force)-2.*stretching_energy(K,l,p,force,fmin)
print(k0*np.exp(-DeltaG/2.*K)*dt)
print(k0*np.exp(DeltaG/2.*K)*dt)
print(np.exp(-DeltaG/2.*K)/(np.exp(DeltaG/2.*K)))
print(np.exp(-DeltaG))
# Calculate random traces with same transition rates 
traces = pd.DataFrame()
cont=0
for m in range(N_traces):
    
    N = np.zeros(N_steps)
    
    for n in range (1,N_steps):
       #probabilities of going forward & backward
        #half the times we check to go forward or backwards
        alea=np.random.random_sample(None)
        alea2=np.random.random_sample(None)
        step=0
        if alea<0.5: #backwards = rezipping
            DeltaG=n_atp*Deltamu-Deltabp+2.*force*xdefWLC(K,l,p,force)-2.*stretching_energy(K,l,p,force,fmin)    
            pre=k0*np.exp(DeltaG/(2.*K))*dt
            #print (pre)
            if alea2>=pre :
                step=-1
                #print (step)
                print ("rezipping")
            else :
                print ("pause")
        if alea>=0.5: #forward = unzipping
            DeltaG=n_atp*Deltamu-Deltabp+2.*force*xdefWLC(K,l,p,force)-2.*stretching_energy(K,l,p,force,fmin)   
            pun=k0*np.exp(-DeltaG/(2.*K))*dt
            #print (pun)
            if alea2>=pun :
                step=1
                cont+=1
                #print (step)
                print ("unzipping")
            else :
                print ("pause")
        #returns 
    #np.random.random(1) returns a random number betweem 0 and 1
    #np.sign returns +/1 depending on the sign inside it
        #step = np.sign(np.random.random(1) - P_b)
        
        N[n] = N[n-1] + step
    
    traces = traces.append(pd.Series(N),ignore_index = True)
    

traces = np.transpose(traces)
print (cont)
dN = [None]*(N_traces)
dN_dt = pd.DataFrame()
#f & gamma changing
for k in range(1,6):
    #N_traces for each force
    for m in range(N_traces):

        n = 0
        dt_1 = 0
        dt_2 = k*5
        dN_temp = np.zeros(N_steps-dt_2)  
        
        while dt_2 < N_steps:
            
            dN_temp[n] = traces[m][int(dt_2)] - traces[m][int(dt_1)]
            dt_1 += 1
            dt_2 += 1
            n += 1
            
        dN[m] = dN_temp
    
    dN_dt = pd.concat([dN_dt,pd.Series([val for sublist in dN for val in sublist])], axis = 1)

# dN_dt is now a matrix of the distributions of probability to move dx at dt, 
# we now need to fit these to a gaussian or just use the histogram data to calculate the large deviation fuction
    
dN_dt.columns = [5*n for n in range(1,6)]

hist, bin_edges = np.histogram(dN_dt[10].dropna(),bins=50)
bin_width = bin_edges[1]-bin_edges[0]
bin_center = bin_edges + bin_width/2
bin_center = bin_center[0:-1]

#%% These are the formulas for the energetics
    
mu = 20
ATP_per_bp = 1
G0 = 2

S_stretch = force/K - (l**2/(2*l-gamma) + 0.25*(gamma-2*l) + gamma**2/(4*l))/(gamma*p)
S_bp = ATP_per_bp*mu - G0

dN = np.diff(N)
dS = dN*S_bp + (dN*gamma)*S_stretch

S = N*S_bp + (N*gamma)*S_stretch
