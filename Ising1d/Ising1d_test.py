from __future__ import print_function
import json
import numpy as np

#Function t0 make the graph
def set_graph(size,d, pbc = True):
    return dict(Name = 'Hypercube', L = size, Dimension = d, Pbc = pbc)
#Function to specify the machine to be used
def set_machine(machname, alpha):
    return dict(Name = machname, Alpha= alpha)
#Function to specify the sampler to be used
def set_sampler(sampname, replicas = 64):
    if sampname == 'MetropolisLocalPt' or sampname == 'MetropolisExchangePt' or sampname == 'MetropolisHamiltonianPt':
        return dict(Name = sampname, Nreplicas = replicas)
    else: return dict(Name = sampname)

#Functions to built the Hamiltonian
def make_single_site_coupling(h, O, L):
    """ Adds a uniform term h sum_i O_i to the hamiltonian, where O is a local op 
    Returns: site list, op list"""
    sites = [[i] for i in range(L)]
    ops = [(h * O.copy()).tolist() for _ in range(L)]
    return sites, ops
def make_nn_coupling(J, O1, O2, L, bc='periodic'):
    """ make a uniform nearest-neighbor coupling list assuming periodic bcs"""
    sites = [[i, (i+1)%L] for i in range(L)]
    ops = [ (J*np.kron(O1, O2)).tolist() for _ in range(L)]
    return sites,ops
#Makes the operators we need in the Spin basis and n
def sigmax():
    return np.array([[0, 1],[1,0]])
def sigmaz():
    return np.array([[1,0], [0, -1]])
def n():
    return (np.identity(2)+sigmaz())/2
#Makes the Hamiltonian, based on how the ising model is done
def get_hilbert(J, h, L):
 #   if h< 2*J:
 #        return dict(Name= 'Spin', S =0.5, TotalSz=0, Nspins=L)
 #   else: 
 #        return dict(Name= 'Spin', TotalSz=0, Nspins=L)
     return dict(QuantumNumbers = [-1,1], Size = L)
def make_tfi_hamiltonian(h, I, J, L):
    """ returns sites and ops lists"""
    sites, ops = make_single_site_coupling(-h, sigmax(), L)
    sites1, ops1 = make_single_site_coupling(-I, sigmaz(), L)
    sites2, ops2 = make_nn_coupling(-J, sigmaz(), sigmaz(), L)
    sites += sites1 
    sites += sites2
    ops += ops1
    ops += ops2
    return sites, ops
#Making the observables
def get_loc_obs_dict(i, O, name=''):
    sites =[[i]]
    ops = [O.tolist()]
    return dict(ActingOn=sites,Operators=ops,Name=name+str(i))
def get_2site_obs_dict(i,j, O1, O2, name = ''):
    sites = [[i, j]]
    ops = [(np.kron(O1, O2)).tolist()]
    return dict(ActingOn=sites, Operators=ops,Name=name+str(i)+str(j))
###############################################################################
#Making the JSON file
###############################################################################
def make_tfi_pars(h, I, J, L, output = 'netket_output'):
    sites, ops = make_tfi_hamiltonian(h, I, J, L)
    #Making the JSON file
    pars = {}
    pars['Graph']    = set_graph(L, 1)
    pars['Hilbert']  = get_hilbert(J, h, L)
    pars['Hamiltonian'] = dict(Operators = ops, ActingOn = sites)
    pars['Machine']  = dict(Name='RbmSpinSymm', Alpha=2.0, UseVisibleBias = True, UseHiddenBias =True)
    pars['Sampler']  = set_sampler('MetropolisLocalPt', replicas = 32)
    pars['Observables'] = [get_2site_obs_dict(i, (i+1)%L, sigmaz(), sigmaz() , name='zz') for i in range(L)]
    for i in range(L):
        pars['Observables'].append(get_loc_obs_dict(i, sigmaz(), name='z'))
        
    pars['Learning'] = {
        'Method'         : 'Sr',
        'Nsamples'       : 1.0e3,
        'NiterOpt'       : 500,
        'Diagshift'      : 0.1,
        'UseIterative'   : False,
        'SaveEvery'      : 10,
        'OutputFile'     : output,
        'StepperType'    : 'Sgd',
        'LearningRate'   : 0.05,
     }
    return pars
