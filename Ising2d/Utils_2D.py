'''These are a series of tools for 2D ED codes using QuSpin'''
import numpy as np

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

def get_sent2(psi, basis, subsys=None, return_rdm=None):
    """Return the entanglement entropy  S_2 of psi, living in basis <basis>, 
    computed in the reduced subsystem specified by subsys
    subsys = list of site labels [0, 1, ..., k] specifying the subsystem. 
    If subsys=None,  defaults to 0....N/2 -1 which correspongs to the upper half
    of the square lattice: 
        [[0,1,2...
          ... LxLy/2] Lxly/2+1,...]
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

'''This function is for the Ising2D model'''
def get_nn_interactions(J,Lx, Ly):
    '''Gets interactions with neares neighbors'''
    nn = []
    
    for i in range(Lx*Ly):
        for j in range(i+1, Lx*Ly):
            if dist(i,j, Lx, Ly) < 1.4:
                nn.append([-J, i, j])
    return nn
'''This function is for the Rydberg2D model'''

def get_kn_interactions(V,ktrunc, Lx, Ly):
    '''Gets interactions with al most k-apart nearest neighbors'''
    nn = []
    
    for i in range(Lx*Ly):
        for j in range(i+1, Lx*Ly):
            d = dist(i,j, Lx, Ly)
            if d < ktrunc+1e-4:
                nn.append([V/(d)**6, i, j])
    return nn

    