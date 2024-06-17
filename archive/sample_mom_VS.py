# test script to figure out bugs in my code
# Apr 2023

import numpy as np
import moments
from copy import deepcopy
import scipy as sp
from scipy.sparse import linalg
from scipy.sparse import coo_matrix

## borrowed directly from https://bitbucket.org/simongravel/moments/src/main/moments/Jackknife.pyx
def python2round(f):
    if round(f + 1) - round(f) != 1:
        return f + abs(f) / f * 0.5
    return round(f)

def index_bis(i, n):
    return int(min(max(python2round(i * n / float(n+1)), 2), n-2))

def calcJK13(n):
    J = np.zeros((n,n-1))
    for i in range(n):
        ibis = index_bis(i + 1, n) - 1
        J[i, ibis] = -(1.+n) * ((2.+i)*(2.+n)*(-6.-n+(i+1.)*(3.+n))-2.*(4.+n)*(-1.+(i+1.)*(2.+n))*(ibis+1.)+(12.+7.*n+n**2)*(ibis+1.)**2) / (2.+n) / (3.+n) / (4.+n)
        J[i, ibis - 1] = (1.+n) * (4.+(1.+i)**2*(6.+5.*n+n**2)-(i+1.)*(14.+9.*n+n**2)-(4.+n)*(-5.-n+2.*(i+1.)*(2.+n))*(ibis+1.)+(12.+7.*n+n**2)*(ibis+1.)**2) / (2.+n) / (3.+n) / (4.+n) / 2.
        J[i, ibis + 1] = (1.+n) * ((2.+i)*(2.+n)*(-2.+(i+1.)*(3.+n))-(4.+n)*(1.+n+2.*(i+1.)*(2.+n))*(ibis+1.)+(12.+7.*n+n**2)*(ibis+1.)**2) / (2.+n) / (3.+n) / (4.+n) / 2.
    return J

def calcD(d):
    data = []
    row = []
    col = []
    # loop over the fs elements:
    for i in range(d):
        if i > 1:
            data.append((i-1) * (d-i))
            row.append(i)
            col.append(i - 1)
        if i < d - 2:
            data.append((i+1) * (d-i-2))
            col.append(i + 1)
            row.append(i)
        if i > 0 and i < d - 1:
            data.append(-2 * i * (d-i-1))
            row.append(i)
            col.append(i)

    return coo_matrix((data, (row, col)), shape=(d, d), dtype='float').tocsc()

def calcS(d, ljk):
    # Computes the jackknife-transformed selection matrix 1
    # for the addition of a single sample
    # arrays for the creation of the sparse (coo) matrix
    # data will have matrix entry, row + column have coordinates
    data = []
    row = []
    col = []
    # loop over the fs elements:
    for i in range(d):
        i_bis = index_bis(i, d - 1) # This picks the second jackknife index 
        i_ter = index_bis(i + 1, d - 1) # This picks the third jackknife index
        # coefficients of the selection matrix
        g1 = i * (d-i) / np.float64(d)
        g2 = -(i+1) * (d-1-i) / np.float64(d)

        if i < d - 1 and i > 0: # First deal with non-fixed variants
            data += [g1 * ljk[i - 1, i_bis - 1], g1 * ljk[i - 1, i_bis - 2],
                    g1 * ljk[i - 1, i_bis], g2 * ljk[i, i_ter - 1],
                    g2 * ljk[i, i_ter - 2], g2 * ljk[i, i_ter]]
            row += 6 * [i]
            col += [i_bis, i_bis - 1, i_bis + 1,
                    i_ter, i_ter - 1, i_ter + 1]
        
        elif i == 0: # g1=0
            data += [g2 * ljk[i, i_ter - 1],
                     g2 * ljk[i, i_ter - 2], g2 * ljk[i, i_ter]]
            row += 3 * [i]
            col += [i_ter, i_ter - 1, i_ter + 1]
        
        elif i == d - 1: # g2=0
            data += [g1 * ljk[i - 1, i_bis - 1], g1 * ljk[i - 1, i_bis - 2],
                     g1 * ljk[i - 1, i_bis]]
            row += 3 * [i]
            col += [i_bis, i_bis - 1, i_bis + 1]

    return coo_matrix((data, (row, col)), shape=(d, d), dtype='float').tocsc()

## code written by me 
def run_mom_iterate(n, s, Nc, theta):
    # matrix to store expected SFAS (# of gens x # of samples)
    mom = np.zeros((len(Nc)+1,n+1),dtype=np.float32)

    # placeholder for SFS of gen k+1
    momkp1 = np.zeros(n+1,dtype=np.float32)

    # computing an array of points at which pop size changes (going from past to present)
    # in the case with two-epoch model, changepoints = [22000, 2000, 1]
    changepoints = len(Nc) - np.concatenate((np.array([0]),np.where(Nc[:-1] != Nc[1:])[0]+1),axis=0)
    changepoints = np.append(changepoints, 1) 

    # initialize the SFS with mutational input in bin 1 (=n x theta/4N_0)
    mom[len(Nc),1] = n*theta/(4*Nc[0])  # singleton input
    
    # only need to do this once - no dependence on N
    J = calcJK13(n)
    S = 0.5 * s * calcS(n+1, J)

    for i in range(len(changepoints)-1):
        # compute this matrix at every changepoint
        D = 0.25/Nc[-changepoints[i]] * calcD(n+1)

        slv = linalg.factorized(sp.sparse.identity(S.shape[0], dtype="float", format="csc") - 0.5 * (D + S))
        Q = sp.sparse.identity(S.shape[0], dtype="float", format="csc") + 0.5 * (D + S)

        # iterate through each gen between changepoints with same D matrix 
        for gen in np.arange(changepoints[i+1]+1,changepoints[i])[::-1]:
            momkp1 = slv(Q.dot(mom[gen+1,]))
            momkp1[0] = momkp1[n] = 0.0

            mom[gen,] = deepcopy(momkp1)

    return mom[:-1,:]     
