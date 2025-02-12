"""

simulate.py
------------

This code simulates a Markov chain given a grid and transition matrix.

Note 1: Loops, conditional statements and functions that return nothing 
terminate by breaking indentation.

Note 2: Python starts indexing at 0.

Note 3: Visit https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html
for information on syntax.

Note 4: Python is a row major language. 

"""

#%%

from numpy import cumsum,nonzero,zeros
from numpy.random import uniform

#%%

"""

                        Simulate Markov chain.

"""

def simulate(grid,pmat,T):
    """
    
    This code simulates a Markov chain given a grid of points from the 
    discretized process and the associated transition matrix.
    
    Input:
        grid : K x N Grid of discretized points.
        pmat : Transition probability matrix.
        seed : Set seed for rng.
        T    : Number of periods to simulate.
        
    Output:
        y : Simulated series.

    """
    
    #%% Initialize.
     
    cmat = cumsum(pmat,axis=1)                                                  # CDF matrix.
    
    if grid.shape[1]%2 == 0:                                                    # Initial state is the mean.
        state0 = int((grid.shape[1]/2)-1)
    else:
        state0 = int((grid.shape[1]-1)/2)
    
    y = zeros((grid.shape[0],T*2))
    
    for i in range(0,T*2):
        y[:,i] = grid[:,state0]
        state1 = cmat[state0,uniform()<=cmat[state0,:]]
        state0 = nonzero(cmat[state0,:]==state1[0])[0][0]
    
    y = y[:,T:T*2]
    
    #%% Output.
    
    return y