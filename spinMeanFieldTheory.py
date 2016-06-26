from realisticModel import *

# Stores all matrices that diagonalize all of the BVK point Hamiltonians
evecs = array([total_bvk_hamiltonian(k1, k2, k3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1) for k1 in range(3) for k2 in range(3) for k3 in range(3)])


# Map from the sublattice, orbital, element, and spin indices to matrix index
def mapping(d, orb, ele, spin):
    return d + orb * 3 + ele * 15 + spin * 30


# Test the mapping function to make sure it is returning the proper indices
for i in range(2):
    for j in range(2):
        for k in range(5):
            for l in range(3):
                print mapping(l, k, j, i)


# Return an element of the matrix that diagonalizes the Hamiltonian
def m(k_point, d, orb, ele, spin):
    return evecs[k_point][mapping(d, orb, ele, spin)]


# Defines a ferromagnetic spin configuration given an Sz value
def sz(sz_val):
    return array([[0, 0, sz_val] for _ in range(27)])


# Defines a ferromagnetic mean field parameter to be used in the self-consistent mean-field theory calculation
def sz_mf(k_point):
    return conj(m(k_point, 0)) * m(k_point, 0) - conj(m(k_point, 1)) * m(k_point, 1)


# Constructs a dynamic spin structure term (next version of code.. need faster language?)
def spin_structure_factor():
    # return m()*m()*m()*m()
    pass
