<<<<<<< HEAD
"""
Created on Mon Jul 30 23:14:02 2018
@author: oscar
"""
'''This code is to build a Ryd1d Hamiltonian with k-th degree interaction'''
import numpy as np
import json

'''Auxiliary Functions'''
def distance(i,j, L):
    dx = abs(i-j)
    dx = np.abs(dx % L)
    return min(dx, L-dx)
def Vlist(V, L, ktrunc):
    l = []
    for i in range(L):
        for j in range(i+1,L):
            d = distance(i,j,L)
            if d<= ktrunc : l.append([V/(d**6), i, j])
    return l
def sigmaz():
    return np.array([[1,0],[0,-1]])
def sigmax():
    return np.array([[0,1], [1,0]])
def n():
    return np.array([[1,0],[0,0]])
def single_site_coupling(O, h, L):
    ops = [(h *O).tolist() for i in range(L)]
    sts = [[i] for i in range(L)]
    return ops, sts
def two_site_interact(O, P, Vlist):
    nn = np.kron(O,P)
    ops = [(Vlist[i][0]*nn).tolist() for i in range(len(Vlist))]
    sts = [[Vlist[i][1], Vlist[i][2]] for i in range(len(Vlist))]
    return ops, sts
def get_loc_obs_dict(i, O, name=''):
    sites =[[i]]
    ops = [O.tolist()]
    return dict(ActingOn=sites,Operators=ops,Name=name+str(i))
def get_2site_obs_dict(i,j, O1, O2, name = ''):
    sites = [[i, j]]
    ops = [(np.kron(O1, O2)).tolist()]
    return dict(ActingOn=sites, Operators=ops,Name=name+str(i)+str(j))

'''Functions to Construct the JSON fIle'''
def get_graph(L):
    return dict(Name = 'Hypercube', L=L, Dimension = 1, Pbc= True)
def get_Hilbert(L):
    return dict(QuantumNumbers = [-1,1], Size = L)
def get_Hamiltonian(L, V, D, O, ktrunc):
    Vs = Vlist(V/4, L, ktrunc)
    Vbar = np.sum([4*Vs[i][0] for i in range(len(Vs))])
    h_x = O/2
    h_z = (D-Vbar)/2
    
    ops, sts = single_site_coupling(-h_x, sigmax(), L)
    ops1, sts1 = single_site_coupling(-h_z, sigmaz(), L)
    ops += ops1
    sts += sts1
    ops2, sts2 = two_site_interact(sigmaz(), sigmaz(), Vs)
    ops += ops2
    sts += sts2
    
    return dict(Operators = ops, ActingOn = sts)
def get_machine(m, alpha):
    return dict(Name = m, Alpha = alpha,InitRandom= False, InitFile = 'InitialParameters.wf')
def get_sampler(Nrep):
    if Nrep == 0: return dict(Name = 'MetropolisHamiltonian')
    else: dict(Name = 'MetropolisHamiltonianPt', Nreplicas=Nrep)
def get_learning(Niter, eta, outfile):
    pars = {'Method' : 'Sr', 
     'Nsamples' : 1.0e3,
     'NiterOpt' : Niter,
     'Diagshift': 0.1,
     'UseIterative' : False,
     'OutputFile'   : outfile,
     'StepperType'  : 'Sgd',
     'LearningRate' : eta,}
    return pars
def get_observables(L, localMag, Corr):
    obs = []
    if localMag == True:
        for i in range(L): obs.append(get_loc_obs_dict(i, sigmaz(), name='z'))
    if Corr == True:
        for i in range(L): obs.append(get_2site_obs_dict(i,(i+1)%L, sigmaz(), sigmaz(), name = 'zz'))
    return obs

def make_Ryd_pars(L, V, D, O, ktrunc, alpha, outfile, m= 'RbmSpinSymm', Nrep = 0, Niter=500, eta = 0.05):
    if ktrunc >L: 
        ktrunc = L
        print('ktrunc larger than L, reshaped to ktrunc=L')
    NetKet = {}
    
    NetKet['Graph']    = get_graph(L)
    NetKet['Hilbert']  = get_Hilbert(L)
    NetKet['Hamiltonian']  = get_Hamiltonian(L,V,D,O, ktrunc) 
    NetKet['Machine']  = get_machine(m, alpha)
    NetKet['Sampler']  = get_sampler(Nrep)
    NetKet['Learning'] = get_learning(Niter, eta, outfile)
    NetKet['Observables'] = get_observables(L, localMag = True, Corr = True)
          
=======
"""
Created on Mon Jul 30 23:14:02 2018
@author: oscar
"""
'''This code is to build a Ryd1d Hamiltonian with k-th degree interaction'''
import numpy as np
import json

'''Auxiliary Functions'''
def distance(i,j, L):
    dx = abs(i-j)
    dx = np.abs(dx % L)
    return min(dx, L-dx)
def Vlist(V, L, ktrunc):
    l = []
    for i in range(L):
        for j in range(i+1,L):
            d = distance(i,j,L)
            if d<= ktrunc : l.append([V/(d**6), i, j])
    return l
def sigmaz():
    return np.array([[1,0],[0,-1]])
def sigmax():
    return np.array([[0,1], [1,0]])
def n():
    return np.array([[1,0],[0,0]])
def single_site_coupling(O, h, L):
    ops = [(h *O).tolist() for i in range(L)]
    sts = [[i] for i in range(L)]
    return ops, sts
def two_site_interact(O, P, Vlist):
    nn = np.kron(O,P)
    ops = [(Vlist[i][0]*nn).tolist() for i in range(len(Vlist))]
    sts = [[Vlist[i][1], Vlist[i][2]] for i in range(len(Vlist))]
    return ops, sts
def get_loc_obs_dict(i, O, name=''):
    sites =[[i]]
    ops = [O.tolist()]
    return dict(ActingOn=sites,Operators=ops,Name=name+str(i))
def get_2site_obs_dict(i,j, O1, O2, name = ''):
    sites = [[i, j]]
    ops = [(np.kron(O1, O2)).tolist()]
    return dict(ActingOn=sites, Operators=ops,Name=name+str(i)+str(j))

'''Functions to Construct the JSON fIle'''
def get_graph(L):
    return dict(Name = 'Hypercube', L=L, Dimension = 1, Pbc= True)
def get_Hilbert(L):
    return dict(QuantumNumbers = [-1,1], Size = L)
def get_Hamiltonian(L, V, D, O, ktrunc):
    Vs = Vlist(V/4, L, ktrunc)
    Vbar = np.sum([4*Vs[i][0] for i in range(len(Vs))])
    h_x = O/2
    h_z = (D-Vbar)/2
    
    ops, sts = single_site_coupling(-h_x, sigmax(), L)
    ops1, sts1 = single_site_coupling(-h_z, sigmaz(), L)
    ops += ops1
    sts += sts1
    ops2, sts2 = two_site_interact(sigmaz(), sigmaz(), Vs)
    ops += ops2
    sts += sts2
    
    return dict(Operators = ops, ActingOn = sts)
def get_machine(m, alpha):
    return dict(Name = m, Alpha = alpha,InitRandom= False, InitFile = 'InitialParameters.wf')
def get_sampler(Nrep):
    if Nrep == 0: return dict(Name = 'MetropolisHamiltonian')
    else: dict(Name = 'MetropolisHamiltonianPt', Nreplicas=Nrep)
def get_learning(Niter, eta, outfile):
    pars = {'Method' : 'Sr', 
     'Nsamples' : 1.0e3,
     'NiterOpt' : Niter,
     'Diagshift': 0.1,
     'UseIterative' : False,
     'OutputFile'   : outfile,
     'StepperType'  : 'Sgd',
     'LearningRate' : eta,}
    return pars
def get_observables(L, localMag, Corr):
    obs = []
    if localMag == True:
        for i in range(L): obs.append(get_loc_obs_dict(i, sigmaz(), name='z'))
    if Corr == True:
        for i in range(L): obs.append(get_2site_obs_dict(i,(i+1)%L, sigmaz(), sigmaz(), name = 'zz'))
    return obs

def make_Ryd_pars(L, V, D, O, ktrunc, alpha, outfile, m= 'RbmSpinSymm', Nrep = 0, Niter=500, eta = 0.05):
    if ktrunc >L: 
        ktrunc = L
        print('ktrunc larger than L, reshaped to ktrunc=L')
    NetKet = {}
    
    NetKet['Graph']    = get_graph(L)
    NetKet['Hilbert']  = get_Hilbert(L)
    NetKet['Hamiltonian']  = get_Hamiltonian(L,V,D,O, ktrunc) 
    NetKet['Machine']  = get_machine(m, alpha)
    NetKet['Sampler']  = get_sampler(Nrep)
    NetKet['Learning'] = get_learning(Niter, eta, outfile)
    NetKet['Observables'] = get_observables(L, localMag = True, Corr = True)
          
>>>>>>> 2dafb372e09419a2337bfb39e5fa60789c50ad21
    return NetKet