'''Implementation of Rydberg 1d model using bosons'''
from quspin.operators import hamiltonian
from quspin.basis import boson_basis_1d
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

def distance(i,j,L):
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

def get_sent(psi, basis, subsys=None, return_rdm=None):
    sdict= basis.ent_entropy(psi, sub_sys_A=subsys,return_rdm=return_rdm, alpha=1.0)
    SA= sdict['Sent_A'] * len(subsys) 
    if return_rdm is not None:
        sdict['Sent_A']=SA        
        return sdict
    return SA

def Ryd_Hamiltonian_bosons(L, Vlist, Delta, Omega, basis):    # PBC
    LF=[[-Delta,i] for i in range(L)]
    TF=[[-Omega/2,i] for i in range(L)]
    static = [['nn', Vlist], ['n', LF], ['+', TF], ['-', TF]]
    dynamic = []
    H=hamiltonian(static,dynamic,basis=basis,check_symm = True, dtype=np.float64)
    energy, psi = H.eigsh(k=1, which='SA')  
    return energy, psi  

#Setting parameters of our Hamiltonian
L = 4
V = 1
ktrunc = 4
D = 0
Vs = Vlist(V,L,ktrunc)
Vbar = np.sum([Vs[i][0] for i in range(len(Vs))])
basis = boson_basis_1d(L, sps =2)

#Moving parameters
O_min, O_max = 0, 1
Nsamples = 10
O = np.linspace(O_min, O_max, Nsamples)

#Defining a bsis and list to store the values of the obsevables' exp vals
GE = []
ent_entropy = []
n1_exp = []
nn12_exp = []
n1_op = n_op(0,basis)
nn12_op = nn_op(0,1, basis)

#Defining a subsystem for the entanglement entropy and the operators
half=[i for i in range(L//2)]
qrtr=[i for i in range(L//4)]

#Mining the information for the cut
for i in range(Nsamples):
    energy0, psi0 = Ryd_Hamiltonian_bosons(L, Vs, D, O[i], basis)  
    GE.append(energy0)
    ent_entropy.append(get_sent(psi0, basis, half))
    n1_exp.append(n1_op.expt_value(psi0))
    nn12_exp.append(nn12_op.expt_value(psi0))

plt.figure()
plt.plot(O, GE)
plt.show()
