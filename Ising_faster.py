""" From "A SURVEY OF COMPUTATIONAL PHYSICS", Python eBook Version
   by RH Landau, MJ Paez, and CC Bordeianu
   Copyright Princeton University Press, Princeton, 2011; Book  Copyright R Landau, 
   Oregon State Unv, MJ Paez, Univ Antioquia, C Bordeianu, Univ Bucharest, 2011.
   Support by National Science Foundation , Oregon State Univ, Microsoft Corp
   
   Adapted by Lev Kaplan 2019"""  

from pylab import *
from random import *

N     = 100                                         # number of spins
J     = 1                        # Exchange energy (J>0 ferromagnetic)
B     = 0.001
k     = 1.                                            # Boltmann constant
Tmax  = 5.                                              # Temperature
Tmin  = 0.1
numTs = 100
Ts = linspace(Tmax,Tmin,numTs)
state = array([+1]*N)                 # spin state: +1 up or -1 (down)  
                                    # start with all spins up            
seed()                                     # Seed random generator

def energy (S):                                 # Method to calc energy
    spinspin = 0.                       # Sum  energy (interaction term only)
    for  i in range(0,N-1):    # i goes from 0 to N-2
        spinspin += S[i]*S[i+1]
    return -J*spinspin - B*sum(S) 
		
ES = energy(state)                  # energy of initial state
mag = sum(state)                 # initial magnetization

avgEls = []
Cls = []
avgMls = [] 
for T in Ts:                                   # Here is the Metropolis algorithm
    teq = 10000 # Number of time-steps to equilibrium
    for j in range (1,t+1):
        r = int(N*random());          # randomly choose which spin to flip
        state[r] *= -1                 # temporarily flip that spin
        #    ET = energy(state)             # finds energy of the test config.
        ET = ES                           # find change in energy
        if (r>=1): ET+=  -2*J*state[r-1]*state[r]
        if (r<=N-2): ET+=  -2*J*state[r]*state[r+1]
        deltamag = 2*state[r]    #find change in magnetization
        
        p = exp((ES-ET)/(k*T))        # test with Boltzmann factor
        if p < random():             # reject change
            state[r] *= -1    # go back and keep previous energy
        else:
            ES = ET       #update energy and magnetization
            mag+=deltamag

    avgE = 0.
    avgE2 = 0.
    avgM = 0.
    t = 10000 # Number of time-steps after equilibrium
    for j in range (1,t+1):
        avgE += ES
        avgE2 += ES**2
        avgM += mag
        
        r = int(N*random());          # randomly choose which spin to flip
        state[r] *= -1                 # temporarily flip that spin
        #    ET = energy(state)             # finds energy of the test config.
        ET = ES                           # find change in energy
        if (r>=1): ET+=  -2*J*state[r-1]*state[r]
        if (r<=N-2): ET+=  -2*J*state[r]*state[r+1]
        deltamag = 2*state[r]    #find change in magnetization
        
        p = exp((ES-ET)/(k*T))        # test with Boltzmann factor
        if p < random():             # reject change
            state[r] *= -1    # go back and keep previous energy
        else:
            ES = ET       #update energy and magnetization
            mag+=deltamag

    avgE /= t
    avgE2 /= t
    avgM /= t

    C = (avgE2 - avgE**2)/(t*T**2)

    avgEls.append(avgE)
    avgMls.append(avgM)
    Cls.append(C)

print('Average energy: ',avgE)
print('Average magnetization: ',avgM)
print('Specific Heat: ',C)

figure()
plot(Ts,avgEls)
title('avgE')
figure()
plot(Ts,avgMls)
title('avgM')
figure()
plot(Ts,Cls)
title('C')