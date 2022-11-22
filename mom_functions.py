#!/usr/bin/env python
# coding: utf-8

# Python file with functions for method of moments framework for getting likelihoods, etc. 

# numerics + rv stuff
import numpy as np
import scipy as sp
from scipy.stats.distributions import chi2
from scipy.sparse import coo_matrix
from scipy.sparse import linalg
from numpy.random import default_rng
# plotting + misc tools
import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches
# import itertools as it
from copy import deepcopy
# import matplotlib.colors as colors
import moments
import warnings
warnings.filterwarnings('error')

# rng setup
rng = default_rng(100496)

# change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["figure.figsize"] = [5, 3.5]
plt.rcParams["figure.dpi"] = 110
plt.rcParams.update({"figure.facecolor": "white"})

# set numpy print option to a more readable format for floats
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

## borrowed directly from https://bitbucket.org/simongravel/moments/src/main/moments/Jackknife.pyx
def python2round(f):
    if round(f + 1) - round(f) != 1:
        return f + abs(f) / f * 0.5
    return round(f)

# The choice i' in n samples that best approximates the frequency of i/(n + 1) is i*n / (n + 1)
def index_bis(i, n):
    return int(min(max(python2round(i * n / float(n+1)), 2), n-2))

# code borrowed from https://bitbucket.org/simongravel/moments/src/main/moments/Jackknife.pyx  
def calcJK13(n):
    J = np.zeros((n,n-1))
    for i in range(n):
        ibis = index_bis(i + 1, n) - 1
        J[i, ibis] = -(1.+n) * ((2.+i)*(2.+n)*(-6.-n+(i+1.)*(3.+n))-2.*(4.+n)*(-1.+(i+1.)*(2.+n))*(ibis+1.)+(12.+7.*n+n**2)*(ibis+1.)**2) / (2.+n) / (3.+n) / (4.+n)
        J[i, ibis - 1] = (1.+n) * (4.+(1.+i)**2*(6.+5.*n+n**2)-(i+1.)*(14.+9.*n+n**2)-(4.+n)*(-5.-n+2.*(i+1.)*(2.+n))*(ibis+1.)+(12.+7.*n+n**2)*(ibis+1.)**2) / (2.+n) / (3.+n) / (4.+n) / 2.
        J[i, ibis + 1] = (1.+n) * ((2.+i)*(2.+n)*(-2.+(i+1.)*(3.+n))-(4.+n)*(1.+n+2.*(i+1.)*(2.+n))*(ibis+1.)+(12.+7.*n+n**2)*(ibis+1.)**2) / (2.+n) / (3.+n) / (4.+n) / 2.
    return J

def calcJK23(n):
    J = np.zeros((n + 1, n - 1))
    for i in range(n + 1):
        ibis = index_bis(i + 1, n) - 1
        if i == n - 1 or i == n:
            ibis = n - 3
        J[i, ibis] = -(1.+n) * ((2.+i)*(2.+n)*(-9.-n+(i+1.)*(3.+n))-2.*(5.+n)*(-2.+(i+1.)*(2.+n))*(ibis+1.)+(20.+9.*n+n**2)*(ibis+1.)**2) / (3.+n) / (4.+n) / (5.+n)
        J[i, ibis - 1] = (1.+n) * (12.+(1.+i)**2*(6.+5.*n+n**2)-(i+1.)*(22.+13.*n+n**2)-(5.+n)*(-8.-n+2.*(i+1.)*(2.+n))*(ibis+1.)+(20.+9.*n+n**2)*(ibis+1.)**2) / (3.+n) / (4.+n) / (5.+n) / 2.
        J[i, ibis + 1] = (1.+n) * ((2.+i)*(2.+n)*(-4.+(i+1.)*(3.+n))-(5.+n)*(n+2.*(i+1.)*(2.+n))*(ibis+1.)+(20.+9.*n+n**2)*(ibis+1.)**2) / (3.+n) / (4.+n) / (5.+n) / 2.
    return J

def calcD(d):
    # res = np.zeros([d, d])
    # # loop over the fs elements:
    # for i in range(d):
    #     if i > 1:
    #         res[i, i - 1] = (i-1) * (d-i)
    #     if i < d - 2:
    #         res[i, i + 1] = (i+1) * (d-i-2)
    #     if i > 0 and i < d - 1:
    #         res[i, i] = -2 * i * (d-i-1)
    # return res
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

# selection with h!=1/2
def calcS2(d, ljk):
    # arrays for the creation of the sparse (coo) matrix
    data = []
    row = []
    col = []
    for i in range(d):
        i_ter = index_bis(i + 1, d - 1)
        i_qua = index_bis(i + 2, d - 1)
        # coefficients
        g1 = (i+1) / np.float64(d) / (d+1.0) * i * (d-i)
        g2 = -(i+1) / np.float64(d) / (d+1.0) * (i+2) * (d-1-i)
        
        if i < d - 1:
            data += [g1 * ljk[i, i_ter - 1], g1 * ljk[i, i_ter - 2],
                     g1 * ljk[i, i_ter], g2 * ljk[i + 1, i_qua - 1],
                     g2 * ljk[i + 1, i_qua - 2], g2 * ljk[i + 1, i_qua]]
            row += 6 * [i]
            col += [i_ter, i_ter - 1, i_ter + 1,
                    i_qua, i_qua - 1, i_qua + 1]
    
        elif i == d - 1: # g2=0
            data += [g1 * ljk[i, i_ter - 1], g1 * ljk[i, i_ter - 2], g1 * ljk[i, i_ter]]
            row += 3 * [i]
            col += [i_ter, i_ter - 1, i_ter + 1]

    return coo_matrix((data, (row, col)), shape=(d, d), dtype='float').tocsc()

def run_mom_iterate_changing(n, s, Nc, mu, misc):
    mom = np.zeros((len(Nc)+1,n+1),dtype=np.float32)
    # momnp1 = np.zeros(n+1)
    momkp1 = np.zeros(n+1,dtype=np.float32)

    changepoints = len(Nc) - np.concatenate((np.array([0]),np.where(Nc[:-1] != Nc[1:])[0]+1),axis=0)
    changepoints = np.append(changepoints, 0)

    mom[len(Nc),1] = n*mu/(4*Nc[0]) # singleton input
    
    # only need to do this once - no dependence on N
    J = calcJK13(n)
    S = 0.5 * s * calcS(n+1, J)

    for i in range(len(changepoints)-1):
        D = 0.25/Nc[len(Nc)-changepoints[i]] * calcD(n+1)

        slv = linalg.factorized(sp.sparse.identity(S.shape[0], dtype="float", format="csc") - 0.5 * (D + S))
        Q = sp.sparse.identity(S.shape[0], dtype="float", format="csc") + 0.5 * (D + S)

        for gen in np.arange(changepoints[i+1],changepoints[i])[::-1]:
            momkp1 = slv(Q.dot(mom[gen+1,]))
            momkp1[0] = momkp1[n] = 0.0

            mom[gen,] = deepcopy(momkp1)

    return mom[:-1,:]           

## creating a log-likelihood function that actually captures the true value using Poisson dist
def get_lp_xl_pois(pxas, g, sXlred, n=2000, cutoff=2):
    """function to compute L(gamma|Xl), where gamma is a range of values and Xl is a given set of freqs"""
    res = np.empty(np.sum((sXlred>=cutoff) & (sXlred<=n-cutoff+1))) #np.empty(len(Xlred))

    # just performing a search in a look-up table
    for idx, i in enumerate(np.where((sXlred>=cutoff) & (sXlred<=n-cutoff+1))[0]):
        res[idx] = -pxas[g][sXlred[i]] + np.log(pxas[g][sXlred[i]])*sXlred[i] - sp.special.gammaln(sXlred[i]+1)
    
    return res
def get_lp_xl_bin(pxas, g, sXlred, n=2000, cutoff=2):
    """function to compute L(gamma|Xl), where gamma is a range of values and Xl is a given set of freqs"""
    res = np.empty(np.sum((sXlred>=cutoff) & (sXlred<=n-cutoff+1))) #np.empty(len(Xlred))

    # just performing a search in a look-up table
    for idx, i in enumerate(np.where((sXlred>=cutoff) & (sXlred<=n-cutoff+1))[0]):
        res[idx] = sp.stats.binom.logpmf(sXlred[i], n, pxas[g][sXlred[i]])
    
    return res

def get_ll_freqdemchanging(s, opts, n=200,):
    selcoef = 2*10**s

    fs = moments.LinearSystem_1D.steady_state_1D(2000, gamma=-selcoef)
    fs = moments.Spectrum(fs)
    fs.integrate(opts['nu'], opts['T'], gamma=-selcoef, dt_fac=0.0005, theta=opts['theta'])
    fs = fs.project([n]) 
    fs[fs<0] = 0
    
    fs = (1 - opts['p_misid']) * fs + opts['p_misid'] * fs[::-1]

    res = (-fs + np.log(fs) * opts['sfs'] - sp.special.gammaln(opts['sfs']+1)).sum()

    return -res

def get_ll_freqagedemchanging(s, opts, n=200):
    selcoef = 10**s

    fsa = run_mom_iterate_changing(n, -2*selcoef/(opts['Nc'][0]/2), opts['Nc']/2, opts['theta'], {})[::-1]
    fsa[fsa<0] = 0

    fsa = (1 - opts['p_misid']) * fsa + opts['p_misid'] * fsa[:,::-1]

    res = np.nansum(-fsa[:-1,1:] + np.log(fsa[:-1,1:]) * opts['sms'][1:,1:] - sp.special.gammaln(opts['sms'][1:,1:]+1))

    return -res


def get_ll_freqconstant(g, opts, n=2000, cutoff=2):
    gamma = 10**g

    fs = moments.LinearSystem_1D.steady_state_1D(2000, gamma=-gamma)
    fs = moments.Spectrum(fs)
    fs.integrate([1], 3, gamma=-gamma, theta=opts['theta']) ## for PReFerSim, we need 0.5Ne instead of Ne
    fs = fs.project([n]) 
    fs[fs<0] = -fs[fs<0]
    
    fs = (1 - opts['p_misid']) * fs + opts['p_misid'] * fs[::-1]

    res = (-fs + np.log(fs) * opts['sfs'] - sp.special.gammaln(opts['sfs']+1)).sum()

    return -res

def get_ll_freqrecconstant(g, opts, n=2000, cutoff=2):
    gamma = 10**g[0]

    fs = moments.LinearSystem_1D.steady_state_1D(2000, gamma=-gamma, h=g[1], theta=opts['theta'],)
    fs = moments.Spectrum(fs)
    fs.integrate([1], 3, gamma=-gamma, theta=opts['theta'], h=g[1]) ## for PReFerSim, we need 0.5Ne instead of Ne
    fs = fs.project([n]) 
    fs[fs<0] = -fs[fs<0]
    
    fs = (1 - opts['p_misid']) * fs + opts['p_misid'] * fs[::-1]

    res = (-fs + np.log(fs) * opts['sfs'] - sp.special.gammaln(opts['sfs']+1)).sum()

    return -res

def get_ll_freqagerecconstant(g, opts, n=2000, cutoff=2):
    gamma = 10**g[0]

    fsa = run_mom_iterate_constantrec(opts['gens'], n, -gamma/opts['N'], opts['N'], opts['theta'], g[1])[::-1]
    fsa[fsa<0] = -fsa[fsa<0]  
    
    fsa = (1 - opts['p_misid']) * fsa + opts['p_misid'] * fsa[::-1]

    res = np.nansum(-fsa[:-1,1:] + np.log(fsa[:-1,1:]) * opts['sms'][1:,1:] - sp.special.gammaln(opts['sms'][1:,1:]+1))
    
    return -res

def get_ll_freqconstantTE(g, opts, n=2000, cutoff=2):
    gamma = 10**g

    fs = moments.LinearSystem_1D.steady_state_1D(5000, gamma=-gamma)
    fs = moments.Spectrum(fs)
    fs.integrate([0.1], 0.25, gamma=-gamma, theta=opts['theta']) 
    fs.integrate([1], 5, gamma=-gamma, theta=opts['theta'])
    fs = fs.project([n]) 
    fs[fs<0] = -fs[fs<0]

    fs = (1 - opts['p_misid']) * fs + opts['p_misid'] * fs[::-1]

    res = (-fs + np.log(fs) * opts['sfs'] - sp.special.gammaln(opts['sfs']+1)).sum()

    return -res

def get_mean_est(g, opts):
    ''' This function was written to check if there was a way to recover the hyperparameter used to simulate
    the expectation (from which we draw the Poisson data). It looks like we can (truth lies within 95% CI). 
    '''
    fsa = rng.poisson(g, size=opts['SMS'].shape)
    fsa[fsa==0] = 1

    res = 0
    for i in range(1,opts['SMS'].shape[0]):
        for j in range(1,opts['SMS'].shape[1]):
            res += -fsa[-i,j] + np.log(fsa[-i,j])*opts['SMS'][i,j] - sp.special.gammaln(opts['SMS'][i,j]+1)

    return -res

def get_ll_freqageconstant(g, opts, n=2000, cutoff=2):
    gamma = 10**g

    fsa = run_mom_iterate_constant(opts['gens'], n, -gamma/opts['N'], opts['N'], opts['theta'], {})[::-1]
    fsa[fsa<0] = -fsa[fsa<0]

    res = np.nansum(-fsa[:-1,1:] + np.log(fsa[:-1,1:]) * opts['sms'][1:,1:] - sp.special.gammaln(opts['sms'][1:,1:]+1))
    
    return -res

def get_ll_freqageconstant_werr(g, opts, n=200,):
    gamma = 10**g

    fsa = run_mom_iterate_constant(opts['gens'], n, -gamma/opts['N'], opts['N'], opts['theta'], {})[::-1]
    fsa[fsa<0] = -fsa[fsa<0]

    # get log-lik from the bins with NO data
    res = -fsa[opts['sms']<1].sum()

    nzidx = opts['sms'].nonzero()
    for i in range(len(nzidx[0])):
        # sd = (np.log(opts['CI'][nzidx[0][i]][nzidx[1][i]][1])-np.log(opts['CI'][nzidx[0][i]][nzidx[1][i]][0]))/(2*1.96)
        # lowci = int(sp.stats.lognorm.ppf(0.01,sd,0,nzidx[0][i])); higci = int(sp.stats.lognorm.ppf(0.99,sd,0,nzidx[0][i]))
        lowci = opts['CI'][nzidx[0][i]][nzidx[1][i]][0]; higci = opts['CI'][nzidx[0][i]][nzidx[1][i]][1]+1 if opts['CI'][nzidx[0][i]][nzidx[1][i]][1]+1 <=opts['gens'] else opts['gens']
        sd = (np.log(higci)-np.log(lowci))/(2*1.96)
        prwts = [sp.stats.lognorm.pdf(x,sd,0,nzidx[0][i]) for x in range(lowci,higci)]
        prwts = prwts/np.sum(prwts)

        tempres = np.zeros(higci-lowci)
        for ia, a in enumerate(range(lowci,higci)):
            tempres[ia] = -fsa[a][nzidx[1][i]] + np.log(fsa[a][nzidx[1][i]])*opts['sms'][a][nzidx[1][i]] - sp.special.gammaln(opts['sms'][a][nzidx[1][i]] + 1) + np.log(prwts[ia])

        res += sp.special.logsumexp(tempres)

    return -res

def get_ll_thetaconstant(g, opts, n=200, cutoff=2):
    """ function to calculate log-lik for single gamma & single theta value with constant pop size """
    gamma, theta = 10**g

    fsa = run_mom_iterate_constant(opts['gens'], n, -gamma/opts['N'], opts['N'], theta, {})[::-1]
    fsa[fsa<0] = -fsa[fsa<0]

    # fsa = (1 - opts['p_misid']) * fsa + opts['p_misid'] * fsa[:,::-1]

    res = np.nansum(-fsa[:-1,1:] + np.log(fsa[:-1,1:]) * opts['sms'][1:,1:] - sp.special.gammaln(opts['sms'][1:,1:]+1))
    
    return -res

def get_ll_theta_bottleneck(g, opts, n=200):
    """ function to calculate log-lik for single gamma and two theta values with two unknown changepoints """
    gamma, theta1, theta2, changept, dur = 10**g

    theta_vec = np.repeat(theta1, opts['sms'].shape[0])
    theta_vec[int(changept):int(changept + dur)] = theta2

    fsa = run_mom_iterate_theta(opts['gens'], n, -gamma/opts['N'], opts['N'], theta_vec, {})[::-1]
    fsa[fsa<0] = -fsa[fsa<0]

    res = np.nansum(-fsa[:-1,1:] + np.log(fsa[:-1,1:]) * opts['sms'][1:,1:] - sp.special.gammaln(opts['sms'][1:,1:]+1))

    return -res

def get_ll_theta_twoepoch(g, opts, n=200):
    """ function to calculate log-lik for single gamma and two theta values with single changepoint """
    gamma, theta1, theta2, changept = 10**g
    # changept = g[-1]

    theta_vec = np.repeat(theta1, opts['sms'].shape[0])
    theta_vec[int(changept):] = theta2

    fsa = run_mom_iterate_theta(opts['gens'], n, -gamma/opts['N'], opts['N'], theta_vec, {})[::-1]
    fsa[fsa<0] = -fsa[fsa<0]

    res = np.nansum(-fsa[:-1,1:] + np.log(fsa[:-1,1:]) * opts['sms'][1:,1:] - sp.special.gammaln(opts['sms'][1:,1:]+1))

    return -res

def get_ll_freqchanging(g, opts, n=1000, cutoff=2):
    alpha, beta = 10**g
    dxs = ((opts['gamma'] - np.concatenate(([opts['gamma'][0]], opts['gamma']))[:-1]) / 2 + (np.concatenate((opts['gamma'], [opts['gamma'][-1]]))[1:] - opts['gamma']) / 2)

    weights = sp.stats.gamma.pdf(opts['gamma'], alpha, scale=beta)
    fs = opts['p_xa_s'][0] * sp.stats.gamma.cdf(opts['gamma'][0],alpha,scale=beta)
    for g, dx, w in zip(opts['gamma'], dxs, weights):
        fs += opts['p_xa_s'][g] * dx * w

    # fs = (1 - opts['p_misid']) * fs + opts['p_misid'] * fs[::-1]

    res = np.nansum(-fs + np.log(fs) * opts['sfs'] - sp.special.gammaln(opts['sfs']+1))

    return -res

def get_ll_freqagechanging(g, opts, n=1000, cutoff=2):
    alpha, beta = 10**g
    dxs = ((opts['gamma'] - np.concatenate(([opts['gamma'][0]], opts['gamma']))[:-1]) / 2 + (np.concatenate((opts['gamma'], [opts['gamma'][-1]]))[1:] - opts['gamma']) / 2)

    weights = sp.stats.gamma.pdf(opts['gamma'], alpha, scale=beta)
    fsa = opts['up_xa_s'][0] * sp.stats.gamma.cdf(opts['gamma'][0],alpha,scale=beta)
    for g, dx, w in zip(opts['gamma'], dxs, weights):
        fsa += opts['up_xa_s'][g] * dx * w
    
    fsa[fsa>1000] = 0

    fsa = (1 - opts['p_misid']) * fsa + opts['p_misid'] * fsa[:,::-1]

    res = np.nansum(-fsa[::-1][:-1,1:] + np.log(fsa[::-1][:-1,1:]) * opts['sms'][1:,1:] - sp.special.gammaln(opts['sms'][1:,1:]+1))

    return -res

## packaging into a function for easy manipulation - iteration implementation 
# input: a (number of gens), n (number of samples), s, N (pop size)
# output: mom (number of sites)
def run_mom_iterate_constant(a, n, s, N, theta, misc):
    mom = np.zeros((a+1,n+1),dtype=np.float32)
    # momnp1 = np.zeros(n+1)
    momkp1 = np.zeros(n+1,dtype=np.float32)

    dt = 1

    D = 0.25/N * calcD(n+1)
    J = calcJK13(n)
    S = 0.5 * s * calcS(n+1, J)

    # if N is same across all gens then only have to do this once
    slv = linalg.factorized(sp.sparse.identity(S.shape[0], dtype="float", format="csc") - dt / 2.0 * (D + S))
    Q = sp.sparse.identity(S.shape[0], dtype="float", format="csc") + dt / 2.0 * (D + S)

    mom[a,1] = n*theta/(4*N) # singleton input

    # going from generation 9 to 0
    for gen in np.arange(a)[::-1]:
        momkp1 = slv(Q.dot(mom[gen+1,]))
        momkp1[0] = momkp1[n] = 0.0

        mom[gen,] = deepcopy(momkp1)

    return mom[:-1,:]           

def run_mom_iterate_constantrec(a, n, s, N, theta, h):
    mom = np.zeros((a+1,n+1),dtype=np.float32)
    # momnp1 = np.zeros(n+1)
    momkp1 = np.zeros(n+1,dtype=np.float32)

    dt = 1

    D = 0.25/N * calcD(n+1)
    J = calcJK13(n)
    J2 = calcJK23(n)
    S = h * s * calcS(n+1, J)
    S2 = s * (1-2*h) * calcS2(n+1, J2)

    # if N is same across all gens then only have to do this once
    slv = linalg.factorized(sp.sparse.identity(S.shape[0], dtype="float", format="csc") - dt / 2.0 * (D + S + S2))
    Q = sp.sparse.identity(S.shape[0], dtype="float", format="csc") + dt / 2.0 * (D + S + S2)

    mom[a,1] = n*theta/(4*N) # singleton input

    # going from generation 9 to 0
    for gen in np.arange(a)[::-1]:
        momkp1 = slv(Q.dot(mom[gen+1,]))
        momkp1[0] = momkp1[n] = 0.0

        mom[gen,] = deepcopy(momkp1)

    return mom[:-1,:]       

def run_mom_iterate_theta(a, n, s, N, theta, misc):
    # assert a == len(theta), "length of mutation rate should be equal to number of gens"
    mom = np.zeros((a+1,n+1),dtype=np.float32)
    # momnp1 = np.zeros(n+1)
    momkp1 = np.zeros(n+1,dtype=np.float32)

    dt = 1

    D = 0.25/N * calcD(n+1)
    J = calcJK13(n)
    S = 0.5 * s * calcS(n+1, J)

    # if N is same across all gens then only have to do this once
    slv = linalg.factorized(sp.sparse.identity(S.shape[0], dtype="float", format="csc") - dt / 2.0 * (D + S))
    Q = sp.sparse.identity(S.shape[0], dtype="float", format="csc") + dt / 2.0 * (D + S)

    ## start with base theta = 1000
    mom[a,1] = n*1000/(4*N) # singleton input

    # going from generation 9 to 0
    for gen in np.arange(a)[::-1]:
        momkp1 = slv(Q.dot(mom[gen+1,]))
        momkp1[0] = momkp1[n] = 0.0

        mom[gen,] = deepcopy(momkp1)

    return mom[:-1,:] * theta[::-1].reshape(-1,1)/1000

## function where each generation was integrated to separately
def run_mom_integrate(a, n, s, N, mu, misc):
    fsmat = np.zeros((a,n+1))
    for idt, dt in enumerate(np.linspace(0.5/N,0.5*a/N,a)[::-1]):
        fs = moments.Spectrum(np.zeros(n + 1))
        fs[1] = 1
        fs.integrate([1], dt, gamma=2*s*N, h=0.5, theta=0, dt_fac=misc['dt_fac'], adapt_dt=misc['adapt_dt'])
        fsmat[idt,:] = n*mu*fs
    return fsmat

## function where each generation is only integrated from previous generation
def run_mom_integrate2(a, n, s, N, mu, misc):
    fsmat = np.zeros((a,n+1))
    dt = 0.5/N
    fs = moments.Spectrum(fsmat[-1,:])
    fs[1] = 1
    fs.integrate([1], dt, gamma=2*s*N, h=0.5, theta=0)
    fsmat[-1,:] = fs
    for idt in np.arange(0,a-1)[::-1]:
        # fs = moments.Spectrum(fsmat[idt+1,:])
        fs.integrate([1], dt, gamma=2*s*N, h=0.5, theta=0, dt_fac=misc['dt_fac'], adapt_dt=misc['adapt_dt'])
        fsmat[idt,:] = fs
    return n*mu*fsmat

## function to obtain the log P(X,|gamma)
def get_lp_xl(p_xa_s, g, sXlred, n=2000, cutoff=20):
    """function to compute L(gamma|Xl), where gamma is a range of values and Xl is a given set of freqs"""
    res = np.empty(np.sum((sXlred>=cutoff) & (sXlred<=n-cutoff+1))) #np.empty(len(Xlred))

    # just performing a search in a look-up table
    for idx, i in enumerate(np.where((sXlred>=cutoff) & (sXlred<=n-cutoff+1))[0]):
        res[idx] = p_xa_s[g][sXlred[i]]
    
    return np.log(res)

def get_lp_xl2(g, sXlred, n=2000, cutoff=20):
    """function to compute L(gamma|Xl), where gamma is a range of values and Xl is a given set of freqs"""
    res = np.empty(np.sum((sXlred>cutoff) & (sXlred<n-cutoff+1))) #np.empty(len(Xlred))

    # ub = np.exp(2.*g)*scipy.special.expi(-2.*g*0.25/N) - scipy.special.expi(2.*g*(1-0.25/N)) - np.exp(2.*g)*(np.log(0.25/N) - np.log(1-0.25/N))
    # lb = np.exp(2.*g)*scipy.special.expi(2.*g*(0.25/N-1)) - scipy.special.expi(2.*g*0.25/N) - np.exp(2.*g)*(np.log(1-0.25/N) - np.log(0.25/N))
    ub = np.exp(2.*g)*sp.special.expi(-2.*g*0.5*cutoff/n) - sp.special.expi(2.*g*(1-0.5*cutoff/n)) - np.exp(2.*g)*(np.log(0.5*cutoff/n) - np.log(1-0.5*cutoff/n))
    lb = np.exp(2.*g)*sp.special.expi(2.*g*(0.5*cutoff/n-1)) - sp.special.expi(2.*g*0.5*cutoff/n) - np.exp(2.*g)*(np.log(1-0.25/n) - np.log(0.5*cutoff/n))
    scalfact = (ub - lb)/np.expm1(2.*g)

    # return a vector...
    for isx, sx in enumerate(np.where((sXlred>cutoff) & (sXlred<n-cutoff+1))[0]):
        res[isx] = (1-np.exp(-2*g*(1-sXlred[sx]/n)))/(sXlred[sx]/n*(1-sXlred[sx]/n)*(1-np.exp(-2*g)))

    return np.log(res/scalfact)

## just doing a lookup of sorts for the right probability
def get_lp_alxl(up_xa_s, g, sXlred, alred, n=2000, cutoff=20):
    # Xsamp = np.arange(1,n)/n
    # sXlred = np.around(Xlred*n).astype(int) # rng.binomial(n, Xlred, len(Xlred))
    res = np.empty(np.sum((sXlred>=cutoff) & (sXlred<=n-cutoff+1)))
    for idx, i in enumerate(np.where((sXlred>=cutoff) & (sXlred<=n-cutoff+1))[0]):
        # if too many gens, then pass in a very low number (like -400.)
        # res[i] = np.log(p_xa_s[g][-int(alred[i]),np.argmin(np.abs(Xlred[i]-Xsamp))+1]) if (int(alred[i]<p_xa_s[g].shape[0])) else -400. 
        try:
            res[idx] = np.log(up_xa_s[g][-int(alred[i]),sXlred[i]]) #if (int(alred[i])<up_xa_s[g].shape[0]) else np.median(np.log(up_xa_s[g][0,:]))
        except RuntimeWarning:
            print(g, sXlred[i], alred[i])
        # if np.isinf(res[idx]):
        #     print(i, Xlred[i], alred[i])

    return res

def get_info_content(p_xa_s, up_xa_s, newdat, num_samps=800, num_sims=16, cutoff=10):
    ci_freq, ci_age = np.zeros(num_sims), np.zeros(num_sims)
    for n in range(num_sims):
        newnewdat = newdat[rng.choice(newdat.shape[0], num_samps, replace=False),:]
        sin_onlyfreq = [np.sum(get_lp_xl(p_xa_s, g1, newnewdat[:,5], cutoff=cutoff)) for g1 in gamma]
        sin_onlyage = [np.sum(get_lp_alxl(up_xa_s, g1, newnewdat[:,5], newnewdat[:,2], cutoff=cutoff)) for g1 in gamma]

        ci_freq[n] = np.abs(get_bfq(sin_onlyfreq-np.max(sin_onlyfreq), gamma)[0])
        ci_age[n] = np.abs(get_bfq(sin_onlyage-np.max(sin_onlyage), gamma)[0])

    return [ci_freq, ci_age]

def get_conf_int(loglik, thresh=2):
    mle = gamma[np.argmax(loglik)]

    if mle==np.min(gamma):
        lower_thresh = gamma[1]  
        return [(np.max(loglik)-thresh-loglik[np.argmax(loglik)+1])/(np.max(loglik) - loglik[np.argmax(loglik)+1])/(mle - lower_thresh)+lower_thresh, mle]
    elif mle==np.max(gamma):
        upper_thresh = gamma[-2] 
        return [mle, -thresh/(loglik[np.argmax(loglik)-1] - np.max(loglik))/(upper_thresh - mle)+mle]
    else:
        lower_thresh = gamma[np.argmax(loglik)+1] 
        upper_thresh = gamma[np.argmax(loglik)-1] 
        return [(np.max(loglik)-thresh-loglik[np.argmax(loglik)+1])*(mle - lower_thresh)/(np.max(loglik) - loglik[np.argmax(loglik)+1])+lower_thresh, -thresh*(upper_thresh - mle)/(loglik[np.argmax(loglik)-1] - np.max(loglik))+mle,]
        # return [(np.max(loglik)-thresh-loglik[np.argmax(loglik)+1])/(np.max(loglik) - loglik[np.argmax(loglik)+1])/(mle - lower_thresh)+lower_thresh, -thresh/(loglik[np.argmax(loglik)-1] - np.max(loglik))/(upper_thresh - mle)+mle]

def get_boot_ci(gamma, p_xa_s, up_xa_s, newdat, nsamps=1000, nboot=20, cutoff=2):
    mle = np.zeros((nboot,2))
    sin_onlyfreq, sin_onlyage = np.zeros(len(gamma)), np.zeros(len(gamma))
    for i in range(nboot):
        newdat = newdat[rng.choice(len(newdat),nsamps,replace=True)]
        for ig, g in enumerate(gamma):
            sin_onlyfreq[ig] = np.sum(get_lp_xl(p_xa_s, g, newdat[:,5], cutoff=cutoff))
            sin_onlyage[ig] = np.sum(get_lp_alxl(up_xa_s, g, newdat[:,5], newdat[:,3], cutoff=cutoff))
        
        mle[i,0] = gamma[np.argmax(sin_onlyfreq)]
        mle[i,1] = gamma[np.argmax(sin_onlyage)]

    return mle            

def get_ci(loglik, gamma, thresh=-1.96):
    ## this function returns the y values for when loglik is -2 (max is 0)
    ## loglik is asymptotically Normal so this is 95% CI
    a, b, c = get_bfq(loglik-np.max(loglik), gamma)
    return [-(-b-np.sqrt(b**2-4*a*(c-thresh)))*0.5/a,-(-b+np.sqrt(b**2-4*a*(c-thresh)))*0.5/a]

def get_bfq(loglik, gamma):
    ## does not work for some reason—wasted multiple hours on it...
    # return np.polynomial.polynomial.Polynomial.fit(-gamma[(ig-3):(ig+3)], loglik[(ig-3):(ig+3)], deg=2)
    igamma = gamma[(np.argmax(loglik)-1):(np.argmax(loglik)+2)] if np.argmax(loglik)>0 else gamma[0:3]
    loglik = loglik[(np.argmax(loglik)-1):(np.argmax(loglik)+2)] if np.argmax(loglik)>0 else loglik[0:3]

    rhs = np.array([np.dot(igamma**2,loglik), np.dot(igamma,loglik), np.sum(loglik)])
    lhs = np.array([[np.sum(igamma**4), np.sum(igamma**3), np.sum(igamma**2)],
    [np.sum(igamma**3), np.sum(igamma**2), np.sum(igamma)],
    [np.sum(igamma**2), np.sum(igamma), len(igamma)]])

    return np.linalg.solve(lhs, rhs) 

# num_sims is number of reps to run to calculate prob
# num_samps is number of rows to resample the big data from
# gamma is np.array of values to calculate over
# thresh is threshold to assign significance
def resample_calculateprob_freq(newdat, gamma, num_sims=16, num_samps=500, thresh=0.05, cutoff=10):
    prob = 0.
    sin_onlyfreq = np.empty(len(gamma))
    dub_onlyfreq = np.zeros((len(gamma),len(gamma)))
    for n in np.arange(num_sims):
        newnewdat = newdat[rng.choice(newdat.shape[0], num_samps, replace=False),:]
        for ig, g in enumerate(gamma):
            # sum log prob for each locus
            sin_onlyfreq[ig] = np.sum(get_lp_xl(g, newnewdat[:,5], cutoff=cutoff))
            for ig2, g2 in enumerate(gamma[0:(ig+1)]):
                dub_onlyfreq[ig, ig2] = np.sum(np.log(0.5*np.exp(get_lp_xl(g, newnewdat[:,5], cutoff=cutoff)) + 0.5*np.exp(get_lp_xl(g2, newnewdat[:,5], cutoff=cutoff))))


        estgonlyfreq = gamma[np.argmax(sin_onlyfreq)]

        # estg1onlyfreq = gamma[np.unravel_index(dub_onlyfreq.argmax(), dub_onlyfreq.shape)[0]]
        # estg2onlyfreq = gamma[np.unravel_index(dub_onlyfreq.argmax(), dub_onlyfreq.shape)[1]]
        estg1onlyfreq, estg2onlyfreq = np.take(gamma, np.unravel_index(np.argmax(np.ma.masked_array(dub_onlyfreq, mask)),dub_onlyfreq.shape))

        lambfreq = 2.*(dub_onlyfreq[gamma==estg1onlyfreq,gamma==estg2onlyfreq] - sin_onlyfreq[gamma==estgonlyfreq])

        if(chi2.sf(lambfreq, 1)<thresh):
            prob += 1.

    return [prob/num_sims, estgonlyfreq, np.array([estg1onlyfreq, estg2onlyfreq])]

# num_samps is number of rows to resample the big data from
# gamma is np.array of values to calculate over
# thresh is threshold to assign significance
def resample_calculateprob_age(newdat, gamma, num_sims=16, num_samps=500, thresh=0.05, cutoff=10):
    prob = 0.

    sin_onlyage = np.empty(len(gamma))
    
    dub_onlyage = np.zeros((len(gamma),len(gamma)))
    for n in np.arange(num_sims):
        newnewdat = newdat[rng.choice(newdat.shape[0], num_samps, replace=False),:]
        for ig, g in enumerate(gamma):
            # sum log prob for each locus
            sin_onlyage[ig] = np.sum(get_lp_alxl(g, newnewdat[:,5], newnewdat[:,2], cutoff=cutoff))
            for ig2, g2 in enumerate(gamma[0:(ig+1)]):
                dub_onlyage[ig, ig2] = np.sum(np.log(0.5*np.exp(get_lp_alxl(g, newnewdat[:,5], newnewdat[:,2], cutoff=cutoff)) + 0.5*np.exp(get_lp_alxl(g2, newnewdat[:,5], newnewdat[:,2], cutoff=cutoff))))

        estgonlyage = gamma[np.argmax(sin_onlyage)]        

        # estg1onlyage = gamma[np.unravel_index(dub_onlyage.argmax(), dub_onlyage.shape)[0]]
        # estg2onlyage = gamma[np.unravel_index(dub_onlyage.argmax(), dub_onlyage.shape)[1]]
        estg1onlyage, estg2onlyage = np.take(gamma, np.unravel_index(np.argmax(np.ma.masked_array(dub_onlyage, mask)),dub_onlyage.shape))

        lambonlyage = 2.*(dub_onlyage[gamma==estg1onlyage,gamma==estg2onlyage] - sin_onlyage[gamma==estgonlyage])

        if(chi2.sf(lambonlyage, 1)<thresh):
            prob += 1.

    return [prob/num_sims, estgonlyage, np.array([estg1onlyage, estg2onlyage])]

def logsumexp(a, axis = None, return_sign = False, keepdims = False):
    """
    log(sum(exp(A))) = M + log(sum(exp(A - M)))
    """

    a_max = np.amax(a, axis=axis, keepdims=True)

    if a_max.ndim > 0:
        a_max[~np.isfinite(a_max)] = 0
    elif not np.isfinite(a_max):
        a_max = 0

    tmp = np.exp(a - a_max)

    # suppress warnings about log of zero
    with np.errstate(divide='ignore'):
        s = np.sum(tmp, axis=axis, keepdims=keepdims)
        if return_sign:
            sgn = np.sign(s)
            s *= sgn  # /= makes more sense but we need zero -> zero
        out = np.log(s)

    if not keepdims:
        a_max = np.squeeze(a_max, axis = axis)
    out += a_max

    if return_sign:
        return out, sgn
    else:
        return out
