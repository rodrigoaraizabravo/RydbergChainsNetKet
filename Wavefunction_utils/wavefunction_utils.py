import numpy as np
import json

'''Here are functions that serve to reconstruct the NetKet wave function and analize it'''
def get_spin(L,i):
    sites = []
    while i > 0:
        if i % 2 == 1: sites.append(1)
        else: sites.append(-1)
        i = int(i/2)
    sites = sites + [-1]*(L-len(sites))
    sites = np.array(sites)[::-1]
    return sites
def get_bin(L,i):
    sites = []
    while i > 0:
        sites.append(i % 2)
        i = int(i/2)
    sites = sites + [0]*(L-len(sites))
    sites = np.array(sites)[::-1]
    return sites
def spin_basis_1d(L):
    '''Construct the spin 1/2 basis with local quantun numbers [-1,1] and chainlength L'''
    basis_dic = {}
    basis = []
    for i in range(2**L):
        sites = get_spin(L,i)
        basis_dic[tuple(sites)] = i
        basis.append(np.asarray(sites))
    return basis_dic, basis
def boson_basis_1d(L):
    '''Construct the boson basis with local quantun numbers [0,1] and chainlength L'''
    basis_dic = {}
    basis = []
    for i in range(2**L):
        sites = get_bin(L,i)
        basis_dic[tuple(sites)] = i
        basis.append(np.asarray(sites))
    return basis_dic, basis
def theta(sigma, b, W):
    '''Computes the thetas for the RBMSpin for a particular sigma'''
    return np.transpose(W).dot(sigma)+b
def theta_symm(S,b,W):
    '''Computes the thetas for the RBMSpin for a particular sigma given
    the matrix S of permutations of a particular spin sigma'''
    th = np.transpose(S.dot(W))
    for col in range(len(b)):
        th[col,:] += b[col]
    return th

def logPsi(sigma, a, b, W):
    '''Reconstructs the log Fourier Coeficcient of RBMSpin for a particular spin sigma'''
    return a.dot(sigma) + np.sum(np.cosh(theta(sigma, b, W)))
def Psi(sigma, a, b, W):
    '''Returns the exponen ot logPsi'''
    return np.exp(logPsi(sigma,a,b,W))

def logPsi_symm(sigma, asymm, bsymm, Wsymm):
    '''Reconstructs the log Fourier Coeficcient of RBMSpinSymm for a particular spin sigma'''
    S = np.array(np.array([np.roll(sigma, i) for i in range(len(sigma))]))
    return asymm*np.sum(sigma) + np.sum(np.cosh(theta_symm(S,bsymm,Wsymm)))
def Psi_symm(sigma,asymm, bsymm, Wsymm):
    '''Returns the exponen ot logPsi_symm'''
    return np.exp(logPsi_symm(sigma,asymm,bsymm,Wsymm))

def load_params(file, which='Last'):
    '''Loads the parameter for an RBMSpin from the .wf file. If which = Last, it will return
    the last RBM, if which = last, it return all the RBMs'''
    with open(file, 'r') as handle:
        json_data = [json.loads(line) for line in handle]
    a_list, b_list, W_list = [], [], []
    
    for data in json_data:
        if data["Machine"]["Name"] != 'RbmSpin': raise 'Expecting RbmSpin Machine' 
        rows = len(data["Machine"]["a"])
        cols = len(data["Machine"]["b"])
        a, b = np.zeros(rows, dtype=complex), np.zeros(cols, dtype=complex)
        W    = np.zeros((rows, cols), dtype=complex)
    
    for i in range(rows):
        xa,ya = data["Machine"]["a"][i][0], data["Machine"]["a"][i][1]
        a[i] = complex(xa,ya)
    for j in range(cols):
        xb,yb = data["Machine"]["b"][j][0], data["Machine"]["b"][j][1]
        b[j] = complex(xb,yb)
    
    for i in range(rows):
        for j in range(cols):
            xw,yw = data["Machine"]["W"][i][j][0], data["Machine"]["W"][i][j][1]
            W[i][j] = complex(xw, yw)
    a_list.append(a)
    b_list.append(b)
    W_list.append(W)
    
    if which == 'Last': return a_list[-1], b_list[-1], W_list[-1]
    if which == 'All': return a_list, b_list, W_list

def load_params_symm(file, which = 'Last'):
    '''Loads the parameter for an RBMSpinSymm from the .wf file. If which = Last, it will return
    the last RBM, if which = last, it return all the RBMs'''
    with open(file, 'r') as handle:
        json_data = [json.loads(line) for line in handle]
    a_list, b_list, W_list = [], [], []
    
    if json_data[0]["Machine"]["Name"]!= 'RbmSpinSymm': 
        raise 'Expecting RbmSpin Symmetric Machine' 
    L = json_data[0]["Machine"]["Nvisible"]
    
    for data in json_data:
        len_b = len(data["Machine"]["bsymm"]) 
        a_symm, b_symm = np.zeros(1, dtype= complex), np.zeros(len_b, dtype=complex)
        W_symm   = np.zeros((L, len_b), dtype=complex)
           
        xa,ya = data["Machine"]["asymm"][0], data["Machine"]["asymm"][1]
        a_symm[0] = complex(xa,ya)
        for j in range(len_b):
            xb,yb = data["Machine"]["bsymm"][j][0], data["Machine"]["bsymm"][j][1]
            b_symm[j] = complex(xb,yb)
            for i in range(L):
                for j in range(len_b):
                    xw,yw = data["Machine"]["Wsymm"][i][j][0], data["Machine"]["Wsymm"][i][j][1]
                    W_symm[i][j] = complex(xw, yw)
        a_list.append(a_symm)
        b_list.append(b_symm)
        W_list.append(W_symm)
    
    if which == 'Last': return a_list[-1], b_list[-1], W_list[-1]
    if which == 'All': return a_list, b_list, W_list
    
def wf(L, a, b, W, basis):
    '''Constructs the wave function for a given w,a,b paramaters and a kind of basis'''
    c = []
    
    for sigma in basis:
        c.append(Psi(sigma, a, b, W))
    #Normalization
    c /= np.sqrt(np.sum(np.abs(np.array(c))**2))
    #Checking
    if np.sum(np.abs(c)**2) - 1 < 1e-3: print("Wave function normalized!")
    else: print("Norm is %f" %np.sum(np.abs(c)**2))
        
    return c

#This function constructs the normalized wave function  
def wf_symm(L, asymm, bsymm, Wsymm, basis):
    '''Constructs the symmetric wave function for a given w,a,b paramaters and a kind of basis'''
    c = []
    
    for sigma in basis:
        c.append(Psi_symm(sigma, asymm, bsymm, Wsymm))
    #Normalization
    c /= np.sqrt(np.sum(np.abs(np.array(c))**2))
    #Checking
    if np.sum(np.abs(c)**2) - 1 < 1e-3: print("Wave function normalized!")
    else: print("Norm is %f" %np.sum(np.abs(c)**2))
    
    return c
'''These functions are auxiliary functions'''
def sigmax():
    return np.array([[0,1],[1,0]])
def sigmaz():
    return np.array([[1,0],[0,1]])
def sigmazsigmaz():
    return np.kron(sigmaz(), sigmaz())
def n():
    return 1/2*(sigmax()+np.identity(2))
def nn():
    return np.kron(n(),n())
def z(i, L):
    op_base = [np.identity(2) for i in range(L)]
    op_base[i] = sigmaz()
    return constr_op(op_base)
def zz(i,j, L):
    op_base = [np.identity(2) for i in range(L)]
    op_base[i], op_base[j] = sigmaz(), sigmaz()
    return constr_op(op_base)
def constr_op(oplist):
    """ given a list of local ops, return the tensor product """
    o = oplist[0]
    for i in range(1, len(oplist)):
        o = np.kron(o, oplist[i])
    return o
def expect_val(op, psi):
    '''Computes the expectation value of an operator on a state psi'''
    bra = np.conj(np.transpose(psi))
    ket = op.dot(psi)
    return np.real(bra.dot(ket))


