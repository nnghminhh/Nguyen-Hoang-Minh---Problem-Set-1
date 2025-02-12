from scipy.stats import norm
from numpy import arange,count_nonzero,expand_dims,identity,linalg,linspace,nonzero,ones,prod,tile,zeros
import numpy as np
import pandas as pd

def rouwenhorst(mu,gamma1,sigma,n):

    """
    
    This function discretizes an AR(1) process.
    
            y(t) = mu + gamma1*y(t-1) + eps(t), eps(t) ~ NID(0,sigma^2)
    
    Input:
        mu    : Intercept of AR(1).
        gamma1  : Persistence of AR(1).
        sigma : Standard deviation of error term.
        n    : Number of states.
        
    Output:
        y    : Grid for the AR(1) process.
        pmat : Transition probability matrix.
        
    """
    #Compute Standard Deviation and Mean:
    ar_sd = np.sqrt((sigma**2)/(1-gamma1**2))
    ar_mean = mu/(1-gamma1)

    #Set transition probabilities
    p = (1 + gamma1) / 2
    q = (1 + gamma1) / 2

    #Define the discrete state space
    y1 = ar_mean - np.sqrt(n-1)*ar_sd
    yn = ar_mean + np.sqrt(n-1)*ar_sd
    y_states = linspace(y1,yn,n)

    #The transition matrix for 2 states
    P_2 = np.array([[p, 1 - p],
                    [1 - q, q]])
    
    pmat = P_2

    for k in range(3, n+1):
        P_k = np.zeros((k, k))
        P_a = np.zeros((k, k))
        P_b = np.zeros((k, k))
        P_c = np.zeros((k, k))
        P_d = np.zeros((k, k))

        P_a[0:k-1, 0:k-1] = p * pmat
        P_b[0:k-1, 1:k] = (1 - p) * pmat
        P_c[1:k, 0:k-1] = (1 - q) * pmat
        P_d[1:k, 1:k] = q * pmat

        P_k = P_a + P_b + P_c + P_d

        P_k[1:-1, :] /= 2

        pmat = P_k

    if count_nonzero(pmat.sum(axis=1)<0.99) > 0:
        raise Exception("Some columns of transition matrix don't sum to 1.") 

    return expand_dims(y_states, axis=0), pmat

mu = 0.5
gamma1 = 0.85
sigma = 1
n = 7

y_states, pmat = rouwenhorst(mu,gamma1,sigma,n)

# Display results
print("Discrete States (y):", y_states)
print("\nTransition Matrix (pmat):\n", pd.DataFrame(pmat, columns=[f"S{i+1}" for i in range(n)], index=[f"S{i+1}" for i in range(n)]))
    
