from wavefunction_utils import spin_basis_1d,load_params_symm,wf_symm
from BruteForce_Entropy import swap_expVal
import numpy as np
L=4
a, b, W = load_params_symm('test2.wf')
basis_dic, basis = spin_basis_1d(L)
psi = wf_symm(L, a, b, W, basis)

ent = -np.log(np.real(swap_expVal(basis_dic, psi, L)))
print('The exact result should be 0.6730 and we get from training %f' %ent)

#testing the classical limit
L=4
psi = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
exVal = np.real(swap_expVal(basis_dic, psi, L))
print('We should get 2, and we get %f' %exVal)
# exact  = 0.679112554431