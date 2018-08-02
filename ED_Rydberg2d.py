from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d
from Utils_2D import get_kn_interactions, get_sent2
import numpy as np 
import json
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

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

def Ryd_Hamiltonian(L, Vlist, h_z, h_x, basis):    # PBC
    LF=[[-h_x,i] for i in range(L)]
    TF=[[-h_z,i] for i in range(L)]
    static = [['zz', Vlist], ['x', LF], ['z', TF]]
    dynamic = []
    H=hamiltonian(static,dynamic,basis=basis,check_symm = True, dtype=np.float64)
    energy, psi = H.eigsh(k=1, which='SA')  
    return energy, psi  

#Setting parameters of our Hamiltonian
Lx = 3
Ly = 3
L = Lx*Ly
V = 1
ktrunc = 10
D = 0
Vs = get_kn_interactions(V, ktrunc, Lx, Ly)
Vbar = np.sum([4*Vs[i][0] for i in range(len(Vs))])
basis = spin_basis_1d(L)

#Moving parameters
O_min, O_max = 0.1, 1
Nsamples = 2
O = np.linspace(O_min, O_max, Nsamples)
h_x = O/2
h_z = (D-Vbar)/2

#Defining a bsis and list to store the values of the obsevables' exp vals
GE = []
ent_entropy = []
z1_exp = []
zz12_exp = []
z1_op = z_op(0,basis)
zz12_op = zz_op(0,1, basis)

#Defining a subsystem for the entanglement entropy and the operators
half=[i for i in range(L//2)]
qrtr=[i for i in range(L//4)]

#Mining the information for the cut
for i in range(Nsamples):
    energy0, psi0 = Ryd_Hamiltonian(L, Vs, h_z, h_x[i], basis)  
    GE.append(energy0)
    ent_entropy.append(get_sent(psi0, basis, half))
    z1_exp.append(z1_op.expt_value(psi0))
    zz12_exp.append(zz12_op.expt_value(psi0))

plt.figure()
plt.plot(O, ent_entropy)
plt.show()
        
