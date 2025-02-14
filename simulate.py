"""
simulate.py
------------
This code simulates a Markov chain given a grid and transition matrix.
Note 1: Loops, conditional statements and functions that return nothing terminate
by breaking indentation.
Note 2: Python starts indexing at 0.
Note 3: Visit https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html for
information on syntax.
Note 4: Python is a row major language.
"""
#%% Imports.
from numpy import cumsum,linspace,nonzero,ones,zeros
from numpy.random import uniform
#%% Simulation function.
def simulate(grid,pmat,T):
    """
    This code simulates a Markov chain given a grid of points from the discretized
    process and the associated transition matrix.
    Input:
    grid : K x N Grid of discretized points.
    pmat : Transition probability matrix.
    seed : Set seed for rng.
    T : Number of periods to simulate.
    Output:
    y : Simulated series.
    """
    #%% Initialize.
    N = grid.shape[1] # Number of states.
    pi0 = cumsum(ones(N)/N) # CDF of uniform distribution for initial state.
    init = linspace(0,N-1,N,endpoint=True) # State indices.
    state0 = int(init[uniform()<=pi0][0]) # Initial state.
    #%% Simulate.
    cmat = cumsum(pmat,axis=1) # CDF matrix.
    y = zeros((grid.shape[0],T)) # Container
    for i in range(0,T): # Simulation.
        y[:,i] = grid[:,state0] # Current state.
        state1 = cmat[state0,uniform()<=cmat[state0,:]] # State next period.
        state0 = nonzero(cmat[state0,:]==state1[0])[0][0] # Update index for next period.
    #y = y[:,T:T] # Burn the first half.
    #%% Output.
    return y