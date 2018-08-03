'''This is a bosonic implementation of the Rydberg model in 2D'''
from quspin.operators import hamiltonian
from quspin.basis import boson_basis_1d
from Utils_2D import get_kn_interactions, get_sent2
import numpy as np 
import json
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def n_op(i,basis,check_symm=True,dtype=np.float64):
    actson = [[1,i]] 
    static = [["n", actson ]]
    dynamic=[]
    return hamiltonian(static,dynamic,basis=basis,dtype=dtype,check_symm=check_symm)

def nn_op(i,j, basis, check_symm=True, dtype=np.float64):
    actson = [[1,i,j]]
    static = [["nn", actson]]
    dynamic = []
    return hamiltonian(static, dynamic, basis=basis, dtype = dtype, check_symm = check_symm)

def Ryd_Hamiltonian_bosons(L, Vlist, Delta, Omega, basis):    # PBC
    LF=[[-Delta,i] for i in range(L)]
    TF=[[-Omega/2,i] for i in range(L)]
    static = [['nn', Vlist], ['n', LF], ['+', TF], ['-', TF]]
    dynamic = []
    H=hamiltonian(static,dynamic,basis=basis,check_symm = True, dtype=np.float64)
    energy, psi = H.eigsh(k=1, which='SA')  
    return energy, psi  

#Setting parameters of our Hamiltonian
Lx = 3
Ly = 3
L = Lx*Ly
V = 1
ktrunc = 3
D = 0
Vs = get_kn_interactions(V, ktrunc, Lx, Ly)
Vbar = np.sum([Vs[i][0] for i in range(len(Vs))])
basis = boson_basis_1d(L, sps =2)
O =0.1

energy0, psi0 = Ryd_Hamiltonian_bosons(L, Vs, D, O, basis)  
print(energy0)

        
