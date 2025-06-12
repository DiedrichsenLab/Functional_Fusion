# Script for importing the MDTB data set from super_cerebellum to general format.
import os
import pandas as pd
import numpy as np
import Functional_Fusion.dataset as ds
import Functional_Fusion.matrix as matrix
import matplotlib.pyplot as plt
import scipy.linalg as sl

## Construct design matrix 
def make_design(hrf,seq,iti):
    h=len(hrf)
    N=len(seq)*iti+len(hrf)-1
    Q=len(seq)
    X=np.zeros((N,Q))
    for i in range(Q):
        start = i*iti
        X[start:start+h,seq[i]]=hrf
    return X

def example_design():
    hrf=[1,0.5]
    X = make_design([1,0.5],[0,1,2,3],1)
    X = sl.block_diag(X,X)
    run = np.array([1,1,1,1,1,2,2,2,2,2])
    X = np.c_[X, matrix.indicator(run)]
    return X 

if __name__ == "__main__":
    X = example_design()
    cond = np.array([1,1,2,2,1,2,2,1])
    C= matrix.indicator(cond,positive=True)
    plt.subplot(2,2,1)
    plt.imshow(X)
    plt.title('Design Matrix')
    plt.subplot(2,2,2)
    plt.imshow(sl.inv(X.T@X))
    plt.title('var(beta)')
    plt.subplot(2,2,3)
    plt.imshow(C)
    plt.title('Constrast')
    Cfull = sl.block_diag(C,np.eye(2))
    Xn = X @ Cfull
    plt.subplot(2,2,4)
    plt.imshow(Xn)

    # Make saome control data 
    Y = np.random.randn(X.shape[0],3)
    beta = sl.pinv(X) @ Y
    data  = beta[:-2,:] # remove intercept regressors  
    nbeta  = ds.optimal_contrast([data], C,X)[0]
    nbeta1 = sl.pinv(Xn)@ Y
    assert np.allclose(nbeta, nbeta1[:-2,:]), "Optimal contrast calculation failed"
    pass 