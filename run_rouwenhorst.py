import gc
gc.set_threshold(1,1,1)

import matplotlib.pyplot as plt
import numpy as np

from rouwenhorst import rouwenhorst
from simulate import simulate

#%% Discretize AR(1).
mu = 0.50000 # Intercept.
rho = 0.85000 # Persistence.
sigma = 1.00000 # Std. dev. of error term.
N = 7 # Number of states.
ar = {'rou': rouwenhorst} # Method.
ar_y = {} # Container.
ar_pi = {} # Container.
# Select method and discretize.
for i in ar:
    ar_y[i], ar_pi[i] = ar[i](mu,rho,sigma,N)

print(ar_y)
print(ar_pi)

#%% Simulate the AR(1)'s based on discretized processes.
seed = 2025 # Seed for random number.
T = 50 # Time horizon.
sim = {}
# Simulate the Markov Chain.
for i in ar_y.keys():
    np.random.seed(seed)
    sim["ar_"+i] = simulate(ar_y[i],ar_pi[i],T)

#%% Plot the time series.
time = range(0,T) # X-axis.
# Plot.
for i in sim:
    plt.figure()
    plt.plot(time,np.squeeze(sim[i]))
    plt.xlabel('Time')
    plt.ylabel('Y')
    plt.title("AR(1): "+i)
    plt.show()
