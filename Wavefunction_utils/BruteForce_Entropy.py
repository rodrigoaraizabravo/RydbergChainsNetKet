'''This code is a brute force implementation of the Swap operator trick to 
compute the entanglement entropy'''

'''We assume that this gets run by a script that has already reconstructed
the entire wave function. We need a reference table'''
import numpy as np
from wavefunction_utils import get_spin

def lookup(s1, s2, s12, s21, basis_dic):
    '''Returns the numbers corresponding to s1, s2, s12, s21'''
    s1_idx = basis_dic[tuple(s1)]
    s2_idx = basis_dic[tuple(s2)]
    s12_idx = basis_dic[tuple(s12)]
    s21_idx = basis_dic[tuple(s21)]
    return s1_idx, s2_idx, s12_idx, s21_idx

def swap_expVal(basis_dic, psi, L):
    '''Expectation value of the swap operator'''
    l = L//2
    e = []
    for i in range(2**l):
        for j in range(2**l):
            for k in range(2**l):
                for p in range(2**l):
                    s1A, s1Ac = get_spin(l,i), get_spin(l,j)
                    s2A, s2Ac = get_spin(l,k), get_spin(l,p)
                    s1 = np.concatenate((s1A, s1Ac))
                    s2 = np.concatenate((s2A, s2Ac))
                    s12 = np.concatenate((s1A, s2Ac))
                    s21 = np.concatenate((s2A, s1Ac))
                    s1_idx, s2_idx, s12_idx, s21_idx = lookup(s1, s2, s12, s21, basis_dic)
                    e.append(psi[s1_idx]*psi[s2_idx]*np.conjugate(psi[s12_idx]*psi[s21_idx]))
    return np.sum(e)


    