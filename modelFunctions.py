from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt
import matplotlib
from cmath import *
import itertools

# define imaginary i
i = complex(0, 1)

# lattice parameters (in Angstroms)
a = 6.0675
c = 3.6570

# primitive vectors
primitive_vectors = array([[a, 0, 0], [-(a/2), a*sqrt(3)/2, 0], [0, 0, c]])

# Deflection angles, I found these angles by hand using VESTA
angle1 = (pi/180) * (180 - 90.0604)
angle2 = (pi/180) * (180 - 75.4757)
angle3 = (pi/180) * (180 - 138.6132)
angle4 = (pi/180) * (180 - 120)
angle5 = (pi/180) * (180 - 100.9823)
angle6 = (pi/180) * (180 - 142.9040)
ang1 = (pi/180) * (180 - 67.7109)
ang2 = (pi/180) * (180 - 70.925)
ang3 = (pi/180) * (180 - 134.9172)
ang4 = (pi/180) * (180 - 71.452)
ang5 = (pi/180) * (180 - 138.6132)


# define reciprocal lattice vectors from the primitive vectors
def reciprocal_vectors(ind1, ind2, ind3):
    cp = cross(primitive_vectors[ind2], primitive_vectors[ind3])
    return 2*pi*cp/dot(primitive_vectors[ind1], cp)


# define a Cartesian k vector to dot with real space vector
def k_vector(kx, ky, kz):
    return array([kx, ky, kz])


# define the fourier phase for reuse in Hamiltonian definitions
def fourier(kx, ky, kz, phase_vector, ang):
    return cos(ang) * exp(i * dot(k_vector(kx, ky, kz), phase_vector))


# Defines the adjoint of a given matrix
def adjoint(c_mat):
    return conj(transpose(c_mat))


# Read data from file
def read_data(data_file):
    f = open(data_file, 'r')
    data = []
    for line in f:
        data.append(line.split())
    return data


# Add the upper triangular half of a matrix to its adjoint to give full Hamiltonian
def full_hamiltonian(ham):
    return ham + adjoint(ham)
