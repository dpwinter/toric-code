import numpy as np

def pdist(x1,x2,L):
    adist = abs(x1-x2)
    return min(adist, L-adist)

def pmanhatten(a,b,L):
    rdist = pdist(a[0],b[0],L)
    cdist = pdist(a[1],b[1],L)
    return rdist + cdist

def pdir(x1,x2,L):
    diff = x2 - x1
    if abs(diff) > L - abs(diff): # via boundary
        return -np.sign(diff)
    else: # via grid
        return np.sign(diff)
