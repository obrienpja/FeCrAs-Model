from numpy import *
from numpy.linalg import *
import matplotlib.pyplot as plt
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
    return c_mat + conj(transpose(c_mat))


# Add the upper triangular half of a matrix to its adjoint to give full Hamiltonian
def full_hamiltonian(ham):
    return ham + adjoint(ham)


# The next 3 definitions construct subsectors of the Cr-Cr sector
def hamiltonian_cr_cr_1(kx, ky, kz, tz):
    ham = zeros((3, 3), dtype='complex')
    ham[0][1] = fourier(kx, ky, kz, -(primitive_vectors[0] + primitive_vectors[1]), angle4)
    ham[0][2] = fourier(kx, ky, kz, -primitive_vectors[0], angle4)
    ham[1][2] = fourier(kx, ky, kz, primitive_vectors[1], angle4)
    ham = -tz * ham
    return ham


def hamiltonian_cr_cr_2(kx, ky, kz, tperp, tzp, t2p):
    ham = zeros((6, 6), dtype='complex')
    # first row
    ham[0][1], ham[0][2] = (cos(angle2) * tperp,) * 2
    ham[0][3] = tzp * fourier(kx, ky, kz, primitive_vectors[2], angle2)
    ham[0][4], ham[0][5] = (t2p * fourier(kx, ky, kz, primitive_vectors[2], angle3),) * 2
    # second row
    ham[1][2] = cos(angle2) * tperp
    ham[1][3], ham[1][5] = (t2p * fourier(kx, ky, kz, primitive_vectors[2], angle3),) * 2
    ham[1][4] = tzp * fourier(kx, ky, kz, primitive_vectors[2], angle1)
    # third row
    ham[2][3], ham[2][4] = (t2p * fourier(kx, ky, kz, primitive_vectors[2], angle3),) * 2
    ham[2][5] = tzp * fourier(kx, ky, kz, primitive_vectors[2], angle1)
    # fourth row
    ham[3][4], ham[3][5] = (cos(angle2) * tperp,) * 2
    # fifth row
    ham[4][5] = cos(angle2) * tperp
    # tack a minus sign on so the hopping processes lower energy
    ham = -ham
    return ham


def hamiltonian_cr_cr_3(kx, ky, kz, tperp, tzp, t2p):
    ham = zeros((6, 6), dtype='complex')
    # first row
    ham[0][1] = tperp * fourier(kx, ky, kz, -primitive_vectors[0], angle2)
    ham[0][2] = tperp * fourier(kx, ky, kz, primitive_vectors[1], angle2)
    ham[0][3] = tzp * fourier(kx, ky, kz, primitive_vectors[2], angle1)
    ham[0][4] = t2p * fourier(kx, ky, kz, -primitive_vectors[0] + primitive_vectors[2], angle3)
    ham[0][5] = t2p * fourier(kx, ky, kz, primitive_vectors[1] + primitive_vectors[2], angle3)
    # second row
    ham[1][2] = tperp * fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1], angle2)
    ham[1][3] = t2p * fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[2], angle3)
    ham[1][4] = tzp * fourier(kx, ky, kz, primitive_vectors[2], angle1)
    ham[1][5] = t2p * fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1] + primitive_vectors[2], angle3)
    # third row
    ham[2][3] = t2p * fourier(kx, ky, kz, -primitive_vectors[1] + primitive_vectors[2], angle3)
    ham[2][4] = t2p * fourier(kx, ky, kz, -primitive_vectors[0] - primitive_vectors[1] + primitive_vectors[2], angle3)
    ham[2][5] = tzp * fourier(kx, ky, kz, primitive_vectors[2], angle1)
    # fourth row
    ham[3][4] = tperp * fourier(kx, ky, kz, -primitive_vectors[0], angle2)
    ham[3][5] = tperp * fourier(kx, ky, kz, primitive_vectors[1], angle2)
    # fifth row
    ham[4][5] = tperp * fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1], angle2)
    # tack a minus sign on so the hopping processes lower energy
    ham = -ham
    return ham


# Constructs the Cr-Cr sector
def hamiltonian_cr_cr(kx, ky, kz, tz, tperp, tzp, t2p):
    ham1 = full_hamiltonian(hamiltonian_cr_cr_1(kx, ky, kz, tz))
    ham2 = full_hamiltonian(hamiltonian_cr_cr_2(kx, ky, kz, tperp, tzp, t2p))
    ham3 = full_hamiltonian(hamiltonian_cr_cr_3(kx, ky, kz, tperp, tzp, t2p))
    slice1 = hstack([ham2, zeros((6, 6), dtype='complex'), zeros((6, 3), dtype='complex')])
    slice2 = hstack([zeros((6, 6), dtype='complex'), ham3, zeros((6, 3), dtype='complex')])
    slice3 = hstack([zeros((3, 6), dtype='complex'), zeros((3, 6), dtype='complex'), ham1])
    return vstack([slice1, slice2, slice3])


# The next 3 definitions construct Fe-Fe subsectors
def hamiltonian_fe_fe_1(kx, ky, kz, tz_fe, tz_fe_p):
    ham = zeros((3, 3), dtype='complex')
    # first row
    ham[0][0] = tz_fe * fourier(kx, ky, kz, primitive_vectors[2], angle5)
    ham[0][1] = tz_fe_p * fourier(kx, ky, kz, primitive_vectors[2], angle6)
    ham[0][2] = tz_fe_p * fourier(kx, ky, kz, -primitive_vectors[0] - primitive_vectors[1] + primitive_vectors[2], angle6)
    # second row
    ham[1][0] = tz_fe_p * fourier(kx, ky, kz, primitive_vectors[2], angle6)
    ham[1][1] = tz_fe * fourier(kx, ky, kz, primitive_vectors[2], angle5)
    ham[1][2] = tz_fe_p * fourier(kx, ky, kz, -primitive_vectors[0] - primitive_vectors[1] + primitive_vectors[2], angle6)
    # third row
    ham[2][0] = tz_fe_p * fourier(kx, ky, kz, -primitive_vectors[0] - primitive_vectors[1] + primitive_vectors[2],
                                  angle6)
    ham[2][1] = tz_fe_p * fourier(kx, ky, kz, -primitive_vectors[0] - primitive_vectors[1] + primitive_vectors[2],
                                  angle6)
    ham[2][2] = tz_fe * fourier(kx, ky, kz, primitive_vectors[2], angle5)
    ham = -ham
    return ham


def hamiltonian_fe_fe_2(kx, ky, kz, tperp_fe):
    ham = zeros((3, 3), dtype='complex')
    ham[0][1] = fourier(kx, ky, kz, -primitive_vectors[0], angle4)
    ham[0][2] = fourier(kx, ky, kz, -primitive_vectors[0], angle4)
    ham[1][2] = fourier(kx, ky, kz, array([0, 0, 0]), angle4)
    ham = -tperp_fe * ham
    return ham


def hamiltonian_fe_fe_3(kx, ky, kz, tperp_fe):
    ham = zeros((3, 3), dtype='complex')
    ham[0][1] = fourier(kx, ky, kz, primitive_vectors[1], angle4)
    ham[0][2] = fourier(kx, ky, kz, array([0, 0, 0]), angle4)
    ham[1][2] = fourier(kx, ky, kz, -primitive_vectors[1], angle4)
    ham = -tperp_fe * ham
    return ham


# Construct the Fe-Fe sector
def hamiltonian_fe_fe(kx, ky, kz, tz_fe, tz_fe_p, tperp_fe):
    sl1 = hstack([zeros((3, 3), dtype='complex'), hamiltonian_fe_fe_1(kx, ky, kz, tz_fe, tz_fe_p), zeros((3, 9), dtype='complex')])
    sl2 = hstack([conj(transpose(hamiltonian_fe_fe_1(kx, ky, kz, tz_fe, tz_fe_p))), zeros((3, 12), dtype='complex')])
    sl3 = hstack([zeros((3, 6), dtype='complex'), hamiltonian_fe_fe_2(kx, ky, kz, tperp_fe), zeros((3, 6), dtype='complex')])
    sl4 = hstack([zeros((3, 9), dtype='complex'), hamiltonian_fe_fe_3(kx, ky, kz, tperp_fe), zeros((3, 3), dtype='complex')])
    sl5 = zeros((3, 15) , dtype='complex')
    ham = vstack([sl1, sl2, sl3, sl4, sl5])
    return ham


# The following 6 definitions construct subsectors of the Cr-Fe sector
def hamiltonian_cr_fe_1(kx, ky, kz, t_cr_fe):
    ham = zeros((3, 3), dtype='complex')
    # first row
    ham[0][0] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[2], ang3)
    ham[0][1] = fourier(kx, ky, kz, primitive_vectors[2], ang2)
    ham[0][2] = fourier(kx, ky, kz, primitive_vectors[2], ang1)
    # second row
    ham[1][0] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[2], ang1)
    ham[1][1] = fourier(kx, ky, kz, primitive_vectors[2], ang3)
    ham[1][2] = fourier(kx, ky, kz, primitive_vectors[2], ang2)
    # third row
    ham[2][0] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[2], ang2)
    ham[2][1] = fourier(kx, ky, kz, primitive_vectors[2], ang1)
    ham[2][2] = fourier(kx, ky, kz, primitive_vectors[2], ang3)
    ham = -t_cr_fe * ham
    return ham


def hamiltonian_cr_fe_2(kx, ky, kz, t_cr_fe):
    ham = zeros((3, 3), dtype='complex')
    # first row
    ham[0][0] = fourier(kx, ky, kz, primitive_vectors[0], ang3)
    ham[0][1] = fourier(kx, ky, kz, array([0, 0, 0]), ang2)
    ham[0][2] = fourier(kx, ky, kz, array([0, 0, 0]), ang1)
    # second row
    ham[1][0] = fourier(kx, ky, kz, primitive_vectors[0], ang1)
    ham[1][1] = fourier(kx, ky, kz, array([0, 0, 0]), ang3)
    ham[1][2] = fourier(kx, ky, kz, array([0, 0, 0]), ang2)
    # third row
    ham[2][0] = fourier(kx, ky, kz, primitive_vectors[0], ang2)
    ham[2][1] = fourier(kx, ky, kz, array([0, 0, 0]), ang1)
    ham[2][2] = fourier(kx, ky, kz, array([0, 0, 0]), ang3)
    ham = -t_cr_fe * ham
    return ham


def hamiltonian_cr_fe_3(kx, ky, kz, t_cr_fe):
    ham = zeros((3, 3), dtype='complex')
    # first row
    ham[0][0] = fourier(kx, ky, kz, primitive_vectors[2], ang3)
    ham[0][1] = fourier(kx, ky, kz, primitive_vectors[1] + primitive_vectors[2], ang2)
    ham[0][2] = fourier(kx, ky, kz, primitive_vectors[2], ang1)
    # second row
    ham[1][0] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[2], ang1)
    ham[1][1] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1] + primitive_vectors[2], ang3)
    ham[1][2] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[2], ang2)
    # third row
    ham[2][0] = fourier(kx, ky, kz, -primitive_vectors[1] + primitive_vectors[2], ang2)
    ham[2][1] = fourier(kx, ky, kz, primitive_vectors[2], ang1)
    ham[2][2] = fourier(kx, ky, kz, -primitive_vectors[1] + primitive_vectors[2], ang3)
    ham = -t_cr_fe * ham
    return ham


def hamiltonian_cr_fe_4(kx, ky, kz, t_cr_fe):
    ham = zeros((3, 3), dtype='complex')
    # first row
    ham[0][0] = fourier(kx, ky, kz, array([0, 0, 0]), ang3)
    ham[0][1] = fourier(kx, ky, kz, primitive_vectors[1], ang2)
    ham[0][2] = fourier(kx, ky, kz, array([0, 0, 0]), ang1)
    # second row
    ham[1][0] = fourier(kx, ky, kz, primitive_vectors[0], ang1)
    ham[1][1] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1], ang3)
    ham[1][2] = fourier(kx, ky, kz, primitive_vectors[0], ang2)
    # third row
    ham[2][0] = fourier(kx, ky, kz, -primitive_vectors[1], ang2)
    ham[2][1] = fourier(kx, ky, kz, array([0, 0, 0]), ang1)
    ham[2][2] = fourier(kx, ky, kz, -primitive_vectors[1], ang3)
    ham = -t_cr_fe * ham
    return ham


def hamiltonian_cr_fe_5(kx, ky, kz, t_cr_fe_p):
    ham = zeros((3, 3), dtype='complex')
    # first row
    ham[0][0] = fourier(kx, ky, kz, array([0, 0, 0]), ang4)
    ham[0][1] = fourier(kx, ky, kz, array([0, 0, 0]), ang4)
    ham[0][2] = fourier(kx, ky, kz, -primitive_vectors[0] - primitive_vectors[1], ang5)
    # second row
    ham[1][0] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1], ang5)
    ham[1][1] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1], ang4)
    ham[1][2] = fourier(kx, ky, kz, array([0, 0, 0]), ang4)
    # third row
    ham[2][0] = fourier(kx, ky, kz, primitive_vectors[0], ang4)
    ham[2][1] = fourier(kx, ky, kz, primitive_vectors[0], ang5)
    ham[2][2] = fourier(kx, ky, kz, -primitive_vectors[1], ang4)
    ham = -t_cr_fe_p * ham
    return ham


def hamiltonian_cr_fe_6(kx, ky, kz, t_cr_fe_p):
    ham = zeros((3, 3), dtype='complex')
    # first row
    ham[0][0] = fourier(kx, ky, kz, primitive_vectors[2], ang4)
    ham[0][1] = fourier(kx, ky, kz, primitive_vectors[2], ang4)
    ham[0][2] = fourier(kx, ky, kz, -primitive_vectors[0] - primitive_vectors[1] + primitive_vectors[2], ang5)
    # second row
    ham[1][0] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1] + primitive_vectors[2], ang5)
    ham[1][1] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1] + primitive_vectors[2], ang4)
    ham[1][2] = fourier(kx, ky, kz, primitive_vectors[2], ang4)
    # third row
    ham[2][0] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[2], ang4)
    ham[2][1] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[2], ang5)
    ham[2][2] = fourier(kx, ky, kz, -primitive_vectors[1] + primitive_vectors[2], ang4)
    ham = -t_cr_fe_p * ham
    return ham


# Constructs the Cr-Fe sector
def hamiltonian_cr_fe(kx, ky, kz, t_cr_fe, t_cr_fe_p):
    sl1 = hstack([zeros((3, 6), dtype='complex'), hamiltonian_cr_fe_1(kx, ky, kz, t_cr_fe), zeros((3, 6), dtype='complex')])
    sl2 = hstack([zeros((3, 6), dtype='complex'), hamiltonian_cr_fe_2(kx, ky, kz, t_cr_fe), zeros((3, 6), dtype='complex')])
    sl3 = hstack([zeros((3, 9), dtype='complex'), hamiltonian_cr_fe_3(kx, ky, kz, t_cr_fe), zeros((3, 3), dtype='complex')])
    sl4 = hstack([zeros((3, 9), dtype='complex'), hamiltonian_cr_fe_4(kx, ky, kz, t_cr_fe), zeros((3, 3), dtype='complex')])
    sl5 = hstack([hamiltonian_cr_fe_5(kx, ky, kz, t_cr_fe_p), hamiltonian_cr_fe_6(kx, ky, kz, t_cr_fe_p), zeros((3, 9), dtype='complex')])
    return vstack([sl1, sl2, sl3, sl4, sl5])


def non_magnetic_hamiltonian(kx, ky, kz, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe):
    sl1 = hstack([hamiltonian_cr_cr(kx, ky, kz, tz, tperp, txp, t2p), hamiltonian_cr_fe(kx, ky, kz, t_cr_fe, t_cr_fe_p)])
    sl2 = hstack([transpose(conj(hamiltonian_cr_fe(kx, ky, kz, t_cr_fe, t_cr_fe_p))), hamiltonian_fe_fe(kx, ky, kz, tz_fe, tz_fe_p, tperp_fe)])
    return vstack([sl1, sl2])


# Import crystal field splitting results (Fix this)
def import_crystal_fields():
    pass

''' to do list
1) construct fe-cr subpieces (Done)
2) construct fe-fe sector (Done)
3) construct fe-cr sector (Done)
4) add cr cfs
5) add fe cfs
6) construct full model with all 1s as hopping parameters (Done)
7) define m() function (Done)
8) define Sz mean field function (Done)
9) write SCMFT algorithm
10) get out hamiltonian notes to remember where things come from
11) think about calculating jij
'''
