from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from Utils_2D import get_nn_interactions, get_sent2
import numpy as np 
import json
import matplotlib.pyplot as plt
import scipy as sc

def z_op(i,basis,check_symm=True,dtype=np.float64):
    actson = [[1,i]] 
    static = [["z", actson ]]
    dynamic=[]
    return hamiltonian(static,dynamic,basis=basis,dtype=dtype,check_symm=check_symm)

def zz_op(i,j, basis, check_symm=True, dtype=np.float64):
    actson = [[1,i,j]]
    static = [["zz", actson]]
    dynamic = []
    return hamiltonian(static, dynamic, basis=basis, dtype = dtype, check_symm = check_symm)

def make_tlfi(L, InterList, I, h, basis):
    h_field=[[-h,i] for i in range(L)]
    I_field=[[-I,i] for i in range(L)]
    static = [['zz', InterList], ['x', h_field], ['z', I_field]]
    dynamic = []
    H=hamiltonian(static,dynamic,basis=basis,check_symm = True, dtype=np.float64)
    energy, psi = H.eigsh(k=1, which='SA')

    return energy[0], psi[:,0]
#------------------------------Creating the Cut--------------------------------
'''The following Code is meant to give us plots for energy, entanglement
   entropy, Z1, Z2, Z1Z1 for a cut in parameter space of the Ising1d model'''
#Defining the parameters for our Hamiltonian
Lx = 3
Ly = 3
L = Lx*Ly
J = 1
I = 0
h = 1.5
basis = spin_basis_1d(L)
z1 = z_op(0,basis)
zz12 = zz_op(0,1, basis)
InterList = get_nn_interactions(J, Lx, Ly)
energy, psi = make_tlfi(L, InterList, I, h, basis)
z1_exp = z1.expt_value(psi)
zz12_exp = zz12.expt_value(psi)
print(energy)
print(z1_exp)
print(zz12_exp)
