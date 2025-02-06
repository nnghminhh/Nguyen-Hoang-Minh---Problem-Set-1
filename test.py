"""

test.py
-------

This code tests all codes.

Note 1: Loops, conditional statements and functions that return nothing 
terminate by breaking indentation.

Note 2: Python starts indexing at 0.

Note 3: Visit https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html
for information on syntax.

Note 4: Python is a row major language. 

"""

#%%

import gc
gc.set_threshold(1,1,1)

import matplotlib.pyplot as plt
import numpy as np

from tauchen import tauchen_ar1
from simulate import simulate

#%%

"""

                            Discretize AR(1).

"""

mu    = 0.00005                                                                 # Intercept.                                                                 # Persistence.
rho = [0.75, 0.85, 0.95, 0.99]
sigma = 1.00000                                                                 # Std. dev. of error term.
N     = 7                                                                       # Number of grid points.
m     = 3                                                                       # For Tauchen (1986).
T = 50
time = range(0,T)


ar    = {'tau': tauchen_ar1} 
ar_y  = {}
ar_pi = {}
def simulate(minh):
    for i in ar:
        if i != "tau":
            ar_y[i], ar_pi[i] = ar[i](mu,minh,sigma,N)
        else:
            ar_y[i], ar_pi[i] = ar[i](mu,minh,sigma,N,m)
            
    del mu,minh,sigma,N,m
    gc.collect()

    seed = 2025


    sim = {}


    for i in ar_y.keys():
        np.random.seed(seed)
        sim["ar_"+i] = simulate(ar_y[i],ar_pi[i],T)
        
    return sim

for minh in rho:
    sim = simulate(minh)
    
    for key in sim:
        plt.plot(time, np.squeeze(sim[key]), label=f"{key}, rho={minh}")

plt.xlabel("Time")
plt.ylabel("Y")
plt.title("AR(1) Simulations for Different Rho Values")
plt.legend()
plt.grid()
plt.show()
