# -*- coding: utf-8 -*-
"""
Created on Wed Jul 25 10:25:22 2018

@author: oscar
"""

'''This code is meant to produce an initialization for our RBM'''
import random
import numpy as np
import json

def RandomGaussian(par, seed, sigma, mu = 0):
    random.seed(seed)
    for i in range(len(par)):
        par[i] = random.gauss(mu, sigma)
    return par

def NoisyW_noBias(nv, alpha, sigma, UseHiddenBias = True, UseVisibleBias = True):
    npar = nv*alpha
    
    if UseVisibleBias == True:
        npar += 1
    if UseHiddenBias  == True:
        npar += alpha
        
    par = np.empty(npar)
    par = RandomGaussian(par, 22818, sigma)
    
    k = 0
    par[k] = 0
    k += 1
    
    for p in range(alpha):
        par[k] = 0
        k += 1
    
    tojson = SetParameters(nv, alpha, UseVisibleBias, UseHiddenBias, par)
    return tojson

def GHZ_noisyBias(nv, alpha, sigma, UseHiddenBias = True, UseVisibleBias = True):
    npar = nv*alpha
    
    if UseVisibleBias == True:
        npar += 1
    if UseHiddenBias  == True:
        npar += alpha
        
    par = np.empty(npar)
    par = RandomGaussian(par, 22818, sigma)
    
    k = 1
    
    for p in range(nv*alpha):
        par[-k] += 5
        k +=1
    
    tojson = SetParameters(nv, alpha, UseVisibleBias, UseHiddenBias, par)
    return tojson

def GHZ_noBias(nv, alpha, mean, sigma, UseHiddenBias = True, UseVisibleBias = True):
    npar = nv*alpha
    
    if UseVisibleBias == True:
        npar += 1
    if UseHiddenBias  == True:
        npar += alpha
        
    par = np.empty(npar)
    par = RandomGaussian(par, 22818, sigma, mu = mean)
    
    k = 0
    par[k] = 0
    k += 1
    
    for p in range(alpha):
        par[k] = 0
        k += 1
    
    for p in range(nv*alpha):
        par[-k] += 0
        k +=1
    
    
    tojson = SetParameters(nv, alpha, UseVisibleBias, UseHiddenBias, par)
    return tojson

def GHZ_LBias(nv, alpha, mean, sigma, UseHiddenBias = True, UseVisibleBias = True):
    npar = nv*alpha
    
    if UseVisibleBias == True:
        npar += 1
    if UseHiddenBias  == True:
        npar += alpha
        
    par = np.empty(npar)
    par = RandomGaussian(par, 22818, sigma, mu = mean)
    
    k = 0
    par[k] -= mean/2
    k += 1
    
    for p in range(alpha):
        par[k] -= mean*nv/2
        k += 1
    
    for p in range(nv*alpha):
        par[k] += mean
        k +=1
    

    tojson = SetParameters(nv, alpha, UseVisibleBias, UseHiddenBias, par)
    return tojson

def SetParameters(nv, alpha, UseVisibleBias, UseHiddenBias, par):
    
    asymm = np.empty(1)
    bsymm = np.empty(alpha)
    Wsymm = np.empty((nv, alpha))
    k=0

    if UseVisibleBias ==  True:
        asymm = par[k]
        k +=1
    else: 
        asymm = np.zeros(1)
        k +=1

    if UseHiddenBias == True:
        for p in range(alpha):
            bsymm[p] = par[k]
            k += 1
    else: 
        bsymm = np.zeros(alpha)
        k += alpha
    
    for i in range(nv):
        for j in range(alpha):
            Wsymm[i,j] = par[k]
            k += 1        
    
    j = {}
    j["Machine"] = {
        "Name"           : "RbmSpinSymm",
        "Nvisible"       : nv,
        "Alpha"          : alpha,
        "UseVisibleBias" : UseVisibleBias,
        "UseHiddenBias"  : UseHiddenBias,
        "asymm"          : asymm.tolist(),
        "bsymm"          : bsymm.tolist(),
        "Wsymm"          : Wsymm.tolist(),
    }
    
    return j

###############################################################################
'''MAIN: Constructs the file to initialize the parameters'''
###############################################################################
def InitParams(L, alpha, sigma, method, Usea = True, Useb = True, mean = 0):

    if method == 'NoisyW_noBias':
        pars = NoisyW_noBias(L, alpha, sigma, Useb, Usea)
    
    if method == 'GHZ_noisyBias':
        pars = GHZ_noisyBias(L, alpha, sigma, Useb, Usea)
    
    if method == 'GHZ_noBias':
        pars = GHZ_noBias(L, alpha, mean, sigma, Useb, Usea)
    
    if method == 'GHZ_LBias':
        pars = GHZ_LBias(L, alpha, mean, sigma, Useb, Usea)
        
    param_file="InitialParameters.wf"
    with open(param_file, 'w') as outfile:
        json.dump(pars, outfile)

    print('Parameters Initialized with method '+ method)
    print('Parameters saved to', param_file)
    
    
    
    
    
    
    
    


    