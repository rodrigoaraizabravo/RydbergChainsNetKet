# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 10:14:36 2018

@author: oscar
"""
import json
from ParameterGenerator import InitParams

def make_par_file(i):
    return 'Ryd1d_{0}.json'.format(i) 
def make_output_file(i):
    return 'Ryd1d_output_{0}'.format(i)
def run_netket(Nsamples):
    import subprocess
    for i in range(Nsamples):
        retcode=subprocess.call(['netket', make_par_file(i)])
L=4
V=1
D = 0.1
Nsamples = 15
O = np.linspace(0,0.5, Nsamples)
ktrunc = 4
alpha = 2

InitParams(L, alpha, 0.1, 'NoisyW_noBias')

def generate_files():
    from Ryd1d_test import make_Ryd_pars
    for i in range(Nsamples):
        pf = make_par_file(i)
        of = make_output_file(i)
        Nk_pars = make_Ryd_pars(L, V, D, O[i], ktrunc, alpha, of)
        with open(pf, 'w') as outfile:
            json.dump(Nk_pars, outfile)

generate_files()
run_netket(Nsamples)
        