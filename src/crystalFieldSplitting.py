from modelFunctions import *


# Construct a unitary matrix to rotate columns of the overlap matrix in Monte Carlo
def unitary(mat_size, m, n, theta):
    uni = identity(mat_size, dtype='complex')
    uni[m][m], uni[n][n] = (cos(theta),)*2
    uni[m][n] = sin(theta)
    uni[n][m] = -sin(theta)
    return uni


# Makes a list of unitary matrices as informed by the changes array (Here I convert from base-1 to base-0 indices)
def unitary_matrices(mat_size, changes_array):
    return [unitary(mat_size, int(changes[0])-1, int(changes[1])-1, float(changes[2])) for changes in changes_array]


# Construct total unitary matrix to change between the Cartesian basis and the Monte Carlo basis
def construct_u_total(mat_size, changes_array):
    u_mat_list = unitary_matrices(mat_size, changes_array)
    u_tot = identity(mat_size)
    for u in u_mat_list:
        u_tot = dot(u, u_tot)
    return u_tot


# This needs to be fixed! Write more functions here!
eigvals_cr =  list(map(float,array(read_data("/home/solidangle/pyramid/CFSResults/eigvalsCr.txt")).flatten()))
eigvals_fe =  list(map(float,array(read_data("/home/solidangle/tetrahedron/CFSResults/eigvalsFe.txt")).flatten()))
eigs_cr = array(read_data("/home/solidangle/pyramid/CFSResults/eigsCr.txt"))
eigs_fe = array(read_data("/home/solidangle/tetrahedron/CFSResults/eigsFe.txt"))
eigenvectors_cr = [map(float, eigs_cr[i]) for i in range(5)]
eigenvectors_fe = [map(float, eigs_fe[i]) for i in range(5)]
changes_cr = read_data("/home/solidangle/pyramid/CFSResults/pyramid_changes_array.txt")
changes_fe = read_data("/home/solidangle/tetrahedron/CFSResults/tetrahedron_changes_array.txt")
pyr_u = construct_u_total(5, changes_cr)
tetra_switch = array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1], [0, 0, 0, 1, 0], [0, 0, 1, 0, 0]])
tetra_u = dot(construct_u_total(5, changes_fe), tetra_switch)
u_cr = dot(eigenvectors_cr, transpose(pyr_u))
u_fe = dot(eigenvectors_fe, transpose(tetra_u))


# Create a matrix with off-diagonal components to move between
def off_diagonal(mat_test, uni):
    return dot(dot(adjoint(uni), diag(mat_test)), uni)


# Import crystal field splitting results (Fix this)
def crystal_field_splitting(split_cr, split_fe):
    split_mat_cr = off_diagonal(split_cr * (eigvals_cr - mean(eigvals_cr)), u_cr)
    u_2_cr = array([[0, 0, 0, 0, 1], [0, 0, 0, 1, 0], [0, 0, 1, 0, 0], [0, 1, 0, 0, 0], [1, 0, 0, 0, 0]])
    cfs_cr = dot(dot(adjoint(u_2_cr), split_mat_cr), u_2_cr)
    split_mat_fe = off_diagonal(split_fe * (eigvals_fe - mean(eigvals_fe)), u_fe)
    u_2_fe = array([[0, 1, 0, 0, 0], [1, 0, 0, 0, 0], [0, 0, 0, 1, 0], [0, 0, 1, 0, 0], [0, 0, 0, 0, 1]])
    cfs_fe = dot(dot(adjoint(u_2_fe), split_mat_fe), u_2_fe)
    return [cfs_cr, cfs_fe]


def full_cfs(cfs_mat):
    return kron(cfs_mat, identity(3))
