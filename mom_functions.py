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
import matplotlib.patches as mpatches
import itertools as it
from copy import deepcopy
import matplotlib.colors as colors
import seaborn
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


# In[4]:


N = 10000
s = -10/N # 25/N -> gamma = 50 - strong selection
mu = 1.25e-8 # human mutation rate
n = 2000 # 2 * # of inds sampled, diploid

# start in generation 10 so generation 11 has all zeros (going back in time)
tot_gen = 10000
time_steps = np.linspace(0, tot_gen-1, 100, dtype=int)

mom = np.zeros((tot_gen+1,n+1))
momnp1 = np.zeros(n+1)
momkp1 = np.zeros((tot_gen+1,n+1))

# double precaution - creating a mask
mk = [False] + [True]*(n-1) + [False]

iter = np.arange(1,n)
iterm1p1 = np.arange(2,n-1)

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

def run_mom_iterate_changing(n, s, Nc, mu, misc):
    mom = np.zeros((len(Nc)+1,n+1))
    # momnp1 = np.zeros(n+1)
    momkp1 = np.zeros(n+1)
    
    changepoints = len(Nc) - np.concatenate((np.array([0]),np.where(Nc[:-1] != Nc[1:])[0]+1),axis=0)
    changepoints = np.append(changepoints, 0)

    mom[len(Nc),1] = 1 # singleton input
    
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

def get_ll_freqchanging(g, opts, sXlred, n=1000, cutoff=2):
    fs = moments.Spectrum(np.zeros(n+1))
    fs[1] = 1
    fs.integrate(opts['nu'], opts['T'], gamma=g*opts['Nc0'], dt_fac=opts['dt_fac'], theta=opts['theta'])
    fs[fs<0] = 1e-250
    pxs = fs/np.sum(fs[np.arange(cutoff,n-cutoff+1)])

    res = np.empty(np.sum((sXlred>cutoff) & (sXlred<n-cutoff+1))) #np.empty(len(Xlred))

    # just performing a search in a look-up table
    for idx, i in enumerate(np.where((sXlred>cutoff) & (sXlred<n-cutoff+1))[0]):
        res[idx] = pxs[sXlred[i]]
    
    return -np.sum(np.log(res))

def get_ll_freqagechanging(g, opts, sXlred, alred, n=1000, cutoff=2):
    pxas = run_mom_iterate_changing(n, 2*g, opts['Nc'][::-1], 1.25e-8, misc = {'dt_fac':0.02, 'adapt_dt':True})
    
    pxas[:,np.arange(cutoff,n-cutoff+1)] = pxas[:,np.arange(cutoff,n-cutoff+1)]/np.sum(pxas[:,np.arange(cutoff,n-cutoff+1)]) 

    res = np.empty(np.sum((sXlred>cutoff) & (sXlred<n-cutoff+1)))
    for idx, i in enumerate(np.where((sXlred>cutoff) & (sXlred<n-cutoff+1))[0]):
        res[idx] = pxas[-int(alred[i]),sXlred[i]]

    return -np.sum(np.log(res))

## packaging into a function for easy manipulation - iteration implementation 
# input: a (number of gens), n (number of samples), s, N (pop size)
# output: mom (number of sites)
def run_mom_iterate_constant(a, n, s, N, mu, misc):
    mom = np.zeros((a+1,n+1))
    # momnp1 = np.zeros(n+1)
    momkp1 = np.zeros(n+1)

    dt = 1

    D = 0.25/N * calcD(n+1)
    J = calcJK13(n)
    S = 0.5 * s * calcS(n+1, J)

    # if N is same across all gens then only have to do this once
    slv = linalg.factorized(sp.sparse.identity(S.shape[0], dtype="float", format="csc") - dt / 2.0 * (D + S))
    Q = sp.sparse.identity(S.shape[0], dtype="float", format="csc") + dt / 2.0 * (D + S)

    iter = np.arange(1,n)
    iterm1p1 = np.arange(2,n-1)

    mom[a,1] = 1 # singleton input

    # going from generation 9 to 0
    for gen in np.arange(a)[::-1]:
        # momkp1[iterm1p1] = 0.25/N * (mom[gen+1,iterm1p1-1] * (iterm1p1-1)*(n-iterm1p1+1) + mom[gen+1,iterm1p1+1] * (iterm1p1+1)*(n-iterm1p1-1) - mom[gen+1,iterm1p1] * 2*iterm1p1*(n-iterm1p1))

        # momkp1[1] = 0.25/N * ((n-2) * 2 * mom[gen+1,2] - 2 * (n-1) * mom[gen+1,1])
        # momkp1[n-1] = 0.25/N * ((n-2) * 2 * mom[gen+1,n-2] - 2 * (n-1) * mom[gen+1,n-1])

        # # notice the difference in indexing for LHS
        # momnp1[np.arange(1,n+1)] = (J @ mom[gen+1,iter])

        # momkp1[iter] += 0.5 * s/(n+1) * (iter * (n+1-iter) * momnp1[iter] - (n-iter) * (iter+1) * momnp1[iter+1])
        
        # momkp1[1:n] = mom[gen+1,1:n] + ((D[1:n,1:n] + S[1:n,1:n]) @ mom[gen+1,1:n])

        momkp1 = slv(Q.dot(mom[gen+1,]))
        momkp1[0] = momkp1[n] = 0.0

        mom[gen,] = deepcopy(momkp1)

    return mom[:-1,:]           

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
            res[idx] = np.log(up_xa_s[g][-int(alred[i]),sXlred[i]]) if (int(alred[i])<up_xa_s[g].shape[0]) else np.median(np.log(up_xa_s[g][0,:]))
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

