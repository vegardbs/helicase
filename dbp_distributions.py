# Experimental parameters

p = 0.75   
K = 4.114
l = 0.66

forces = [ 5. ,6. ,7. , 8., 9., 10., 11., 12.]
gammas = [ 0.5538,    0.6155,    0.6668  ,  0.7100 ,   0.7466   , 0.7781  ,  0.8054 ,  0.8293]   
gammas_phi = [  0.3532,    0.4724,    0.5741,    0.6537,    0.7139,    0.7596,    0.7951,    0.8293]


time = 100e-3 

dt = 1e-7*15000
step_dt = int(np.floor(time/dt))

forces = [5,6,7,8,9,10,11,12]

#%%
slopes = np.zeros(len(forces))
N_atp = np.zeros(len(forces))

for k in range(len(forces)):
    
    data = np.loadtxt('/home/vegard/Helicase/python/Simulation_221018_3/' + str(forces[k]) + '/0.dat')
    data2 = np.asarray(data)
    dx = data2[:,1]
    dt = data2[:,0]*15000
    
    plt.plot(dt,dx)
    
    n = 0
    dt_1 = 0
    dt_2 = step_dt
    N_steps =int(len(dt)-dt_2)
    dbp = np.zeros(N_steps)  
    
    while dt_2 < N_steps:
        
        dbp[n] = dx[int(dt_2)] - dx[int(dt_1)]
        dt_1 += 1
        dt_2 += 1
        n += 1
            
    
    slopes[k] = 2*np.mean(dbp)/np.var(dbp)
    N_atp[k] = (slopes[k] - forces[k]/K + (1/p)*(l**2/(2*l-gammas[k]) - (gammas[k]+2*l)/4 + gammas[k]**2/(4*l)) )/20 + 2/20
