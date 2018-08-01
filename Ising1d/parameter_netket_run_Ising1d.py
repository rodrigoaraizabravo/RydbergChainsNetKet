import numpy as np
import json

def make_param_file_name(i):
#    return 'ising1d_test_{0}.json'.format(i)
     return 'ising1d_AnegI01_{0}.json'.format(i)
def make_output_file_name(i):
#    return 'ising1d_test_output_{0}'.format(i)
     return 'ising1d_AnegI01_output_{0}'.format(i)
#Runs Netket on the files
def run_netket(Nsamples):
    import subprocess
    for i in range(Nsamples):
        retcode=subprocess.call([ 'mpirun', '-np', '8', './netket', make_param_file_name(i)] )
######################################################################
#Let us try this on an 1D Ising transverse plus lonfitudinal field Hamiltonian
L=12
J=-1
I = -0.1
hmin = 0
hmax = 1.5
Nsamples = 30
hsamples = np.linspace(hmin, hmax, Nsamples)

# This function generates the parameter files to be run by Netket
def gen_param_files():
    from ising1d_test import make_tfi_pars
    for i in range(Nsamples):
        pfile = make_param_file_name(i)
        outfile = make_output_file_name(i)
        pars = make_tfi_pars(hsamples[i], I, J, L, outfile)
        with open(pfile, 'w') as pf:
            json.dump(pars, pf)

gen_param_files()
run_netket(Nsamples)
