#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 12:14:52 2018

@author: vegard
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



# Experimental parameters

p = 0.75;   
K = 4.114;
l = 0.66; 
forces = [5,6,7,8,9,10,11,12]
gammas = [ 0.5538,    0.6155,    0.6668  ,  0.7100 ,   0.7466   , 0.7781  ,  0.8054 ,  0.8293];    
gammas_phi = [  0.3532,    0.4724,    0.5741,    0.6537,    0.7139,    0.7596,    0.7951,    0.8293];
 
force = forces[6]
gamma = gammas[6]
   
N_steps = 500
N_traces = 50

# Calculate random traces with same transition rates 
traces = pd.DataFrame()

for m in range(N_traces):
    
    N = np.zeros(N_steps)
    
    for n in range (1,N_steps):
        P_f = 0.8
        P_b = 0.3
        step = np.sign(np.random.random(1) - P_b)
        
        N[n] = N[n-1] + step
    
    traces = traces.append(pd.Series(N),ignore_index = True)
    

traces = np.transpose(traces)

dN = [None]*(N_traces)
dN_dt = pd.DataFrame()

for k in range(1,6):
    
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