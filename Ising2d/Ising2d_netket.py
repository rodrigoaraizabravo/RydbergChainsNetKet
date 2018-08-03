'''This code is to build a Ising Hamiltonian in 2d'''
import numpy as np
import json

def get_coordinates(site, Lx, Ly):
    return np.array([site%Lx, site//Ly])
def dmod(dx,L):
    dx = np.abs(dx % L)
    return min(dx, L-dx)
def dist(s1, s2, Lx, Ly):
    x_y = get_coordinates(s1, Lx, Ly) - get_coordinates(s2, Lx, Ly)
    d_c = np.sqrt( (x_y[0])**2 + (x_y[1])**2 )
    dmod2 = np.sqrt( (dmod(x_y[0], Lx))**2 + (dmod(x_y[1], Ly))**2)
    return min(d_c, dmod2)
def Jlist(J, Lx, Ly):
    l = []
    for i in range(Lx*Ly):
        for j in range(i+1,Lx*Ly):
            d = dist(i,j,Lx, Ly)
            if d < 1.4 : l.append([-J, i, j])
    return l
def sigmaz():
    return np.array([[1,0],[0,-1]])
def sigmax():
    return np.array([[0,1], [1,0]])

def single_site_coupling(O, h, L):
    ops = [(h *O).tolist() for i in range(L)]
    sts = [[i] for i in range(L)]
    return ops, sts
def two_site_interact(O, P, Jlist):
    nn = np.kron(O,P)
    ops = [(Jlist[i][0]*nn).tolist() for i in range(len(Jlist))]
    sts = [[Jlist[i][1], Jlist[i][2]] for i in range(len(Jlist))]
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
def get_Hamiltonian(Lx, Ly, J, h, I):
    Js = Jlist(J, Lx, Ly)
    
    ops, sts = single_site_coupling(-h, sigmax(), Lx*Ly)
    ops1, sts1 = single_site_coupling(-I, sigmaz(), Lx*Ly)
    ops += ops1
    sts += sts1
    ops2, sts2 = two_site_interact(sigmaz(), sigmaz(), Js)
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

def make_i2d_pars(Lx,Ly, J, h, I, alpha, outfile, m= 'RbmSpinSymm', Nrep = 0, Niter=1000, eta = 0.05):
    L = Lx*Ly
    NetKet = {}
    
    NetKet['Graph']    = get_graph(L)
    NetKet['Hilbert']  = get_Hilbert(L)
    NetKet['Hamiltonian']  = get_Hamiltonian(Lx,Ly,J,h,I) 
    NetKet['Machine']  = get_machine(m, alpha)
    NetKet['Sampler']  = get_sampler(Nrep)
    NetKet['Learning'] = get_learning(Niter, eta, outfile)
    NetKet['Observables'] = get_observables(L, localMag = True, Corr = True)
    return NetKet

#TEST CODE
parsfile = 'ising2d.json'
Lx, Ly = 3,3
J= 1
h=1.5
I=0
from ParameterGenerator import InitParams
InitParams(Lx*Ly, 4, 0.1, 'NoisyW_noBias')

pars = make_i2d_pars(Lx, Ly, J, h, I, 4, 'test')
with open(parsfile, 'w') as outfile:
    json.dump(pars, outfile)
    
    