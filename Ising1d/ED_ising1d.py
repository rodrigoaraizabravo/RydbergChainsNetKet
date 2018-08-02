<<<<<<< HEAD
from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d, boson_basis_1d
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

def get_sent2(psi, basis, subsys=None, return_rdm=None):
    """Return the entanglement entropy  S_2 of psi, living in basis <basis>, 
    computed in the reduced subsystem specified by subsys
    subsys = list of site labels [0, 1, ..., k] specifying the subsystem. 
    If subsys=None,  defaults to 0....N/2 -1
    
    return_rdm can be specified as 'A' (the subsystem of interest), 
    'B', or both; if so a dictionary is returned
    """
    if subsys is None:
        #the default quspin block
        subsys=tuple(range(basis.N//2))
    
    sdict= basis.ent_entropy(psi, sub_sys_A=subsys,return_rdm=return_rdm, alpha=2.0)
    # the quspin value is normalized by the subsystem size
    SA= sdict['Sent_A'] * len(subsys) 
    if return_rdm is not None:
        sdict['Sent_A']=SA        
        return sdict
    return SA

def make_tlfi(L, J, I, h, basis):
#    J, I, h  = 4*J, 2*I, 2*h
    #Constructing the Hamiltonian
    J_zz=[[-J,i,(i+1)%L] for i in range(L)] # PBC
    h_field=[[-h,i] for i in range(L)]
    I_field=[[-I,i] for i in range(L)]
    static = [['zz', J_zz], ['x', h_field], ['z', I_field]]
    dynamic = []
    H=hamiltonian(static,dynamic,basis=basis,check_symm = True, dtype=np.float64)
    energy, psi = H.eigsh(k=1, which='SA')

    return energy[0], psi[:,0]

#------------------------------Creating the Cut--------------------------------
'''The following Code is meant to give us plots for energy, entanglement
   entropy, Z1, Z2, Z1Z1 for a cut in parameter space of the Ising1d model'''
#Defining the parameters for our Hamiltonian
L = 20
J = 1  
I = 0
h = 1
basis = spin_basis_1d(L)
energy, psi = make_tlfi(L, J, I, h, basis)
print(energy)
=======
from quspin.operators import hamiltonian
from quspin.basis import spin_basis_1d, boson_basis_1d
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

def get_sent2(psi, basis, subsys=None, return_rdm=None):
    """Return the entanglement entropy  S_2 of psi, living in basis <basis>, 
    computed in the reduced subsystem specified by subsys
    subsys = list of site labels [0, 1, ..., k] specifying the subsystem. 
    If subsys=None,  defaults to 0....N/2 -1
    
    return_rdm can be specified as 'A' (the subsystem of interest), 
    'B', or both; if so a dictionary is returned
    """
    if subsys is None:
        #the default quspin block
        subsys=tuple(range(basis.N//2))
    
    sdict= basis.ent_entropy(psi, sub_sys_A=subsys,return_rdm=return_rdm, alpha=2.0)
    # the quspin value is normalized by the subsystem size
    SA= sdict['Sent_A'] * len(subsys) 
    if return_rdm is not None:
        sdict['Sent_A']=SA        
        return sdict
    return SA

def make_tlfi(L, J, I, h, basis):
#    J, I, h  = 4*J, 2*I, 2*h
    #Constructing the Hamiltonian
    J_zz=[[-J,i,(i+1)%L] for i in range(L)] # PBC
    h_field=[[-h,i] for i in range(L)]
    I_field=[[-I,i] for i in range(L)]
    static = [['zz', J_zz], ['x', h_field], ['z', I_field]]
    dynamic = []
    H=hamiltonian(static,dynamic,basis=basis,check_symm = True, dtype=np.float64)
    energy, psi = H.eigsh(k=1, which='SA')

    return energy[0], psi[:,0]

#------------------------------Creating the Cut--------------------------------
'''The following Code is meant to give us plots for energy, entanglement
   entropy, Z1, Z2, Z1Z1 for a cut in parameter space of the Ising1d model'''
#Defining the parameters for our Hamiltonian
L = 20
J = 1  
I = 0
h = 1
basis = spin_basis_1d(L)
energy, psi = make_tlfi(L, J, I, h, basis)
print(energy)
>>>>>>> 2dafb372e09419a2337bfb39e5fa60789c50ad21
