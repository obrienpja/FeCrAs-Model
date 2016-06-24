from realisticModel import *


# Defines Born Von-Karman point for a given system size
def bvk(n1, n2, n3, d1, d2, d3):
    b1 = reciprocal_vectors(0, 1, 2)
    b2 = reciprocal_vectors(1, 2, 0)
    b3 = reciprocal_vectors(2, 0, 1)
    return b1 * float(n1) / d1 + b2 * float(n2) / d2 + b3 * float(n3) / d3


# Constructs the non-magnetic Hamiltonian
def non_magnetic_bvk_hamiltonian(n1, n2, n3, d1, d2, d3, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe):
    k_ptx, k_pty, k_ptz = bvk(n1, n2, n3, d1, d2, d3)
    ham = kron(identity(2), non_magnetic_hamiltonian(k_ptx, k_pty, k_ptz, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe))
    return ham


# Constructs magnetic term of the Hamiltonian
def magnetic_hamiltonian(mu_val, delta, j_h, sz_cr, sz_fe):
    mu_cr = mu_val + delta
    mu_fe = mu_val - delta
    mat1 = diag([-mu_cr - j_h * sz_cr, -mu_fe - j_h * sz_fe, -mu_cr + j_h * sz_cr, -mu_fe + j_h * sz_fe])
    ham = kron(mat1, identity(15))


# Constructs the total Hamiltonian at a given BVK point
def total_bvk_hamiltonian(n1, n2, n3, d1, d2, d3, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe):
    return non_magnetic_bvk_hamiltonian(n1, n2, n3, d1, d2, d3, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe) + magnetic_hamiltonian(mu_val, delta, j_h, sz_cr, sz_fe)


# Returns the matrix that diagonalizes the Hamiltonian
def bvk_eigenvectors(n1, n2, n3, d1, d2, d3, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe):
    ham = total_bvk_hamiltonian(n1, n2, n3, d1, d2, d3, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe)
    return eig(ham)[1]


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
    return conj(m(k_point, 0)) * m(k_point, 0) - m(k_point, 1) * m(k_point, 1)


# Constructs a dynamic spin structure term (next version of code.. need faster language?)
def spin_structure_factor():
    # return m()*m()*m()*m()
    pass
