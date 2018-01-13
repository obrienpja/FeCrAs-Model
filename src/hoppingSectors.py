from modelFunctions import *

# The next 3 definitions construct subsectors of the Cr-Cr sector
def hamiltonian_cr_cr_1(kx, ky, kz, t_a_cr):
    ham = zeros((3, 3), dtype='complex')
    ham[0][1] = fourier(kx, ky, kz, (-primitive_vectors[0] - primitive_vectors[1]), angle4)
    ham[0][2] = fourier(kx, ky, kz, -primitive_vectors[0], angle4)
    ham[1][2] = fourier(kx, ky, kz, primitive_vectors[1], angle4)
    return t_a_cr * ham


def hamiltonian_cr_cr_2(kx, ky, kz, t_in_cr, t_z_cr, t_out_cr):
    ham = zeros((6, 6), dtype='complex')
    # first row
    ham[0][1], ham[0][2] = (t_in_cr * cos(angle2),) * 2
    ham[0][3] = t_z_cr * fourier(kx, ky, kz, primitive_vectors[2], angle1)
    ham[0][4], ham[0][5] = (t_out_cr * fourier(kx, ky, kz, primitive_vectors[2], angle3),) * 2
    # second row
    ham[1][2] = t_in_cr * cos(angle2)
    ham[1][3]  = t_out_cr * fourier(kx, ky, kz, primitive_vectors[2], angle3)
    ham[1][4] = t_z_cr * fourier(kx, ky, kz, primitive_vectors[2], angle1)
    ham[1][5] = t_out_cr * fourier(kx, ky, kz, primitive_vectors[2], angle3)
    # third row
    ham[2][3] = t_out_cr * fourier(kx, ky, kz, primitive_vectors[2], angle3)
    ham[2][4] = t_out_cr * fourier(kx, ky, kz, primitive_vectors[2], angle3)
    ham[2][5] = t_z_cr * fourier(kx, ky, kz, primitive_vectors[2], angle1)
    # fourth row
    ham[3][4], ham[3][5] = (t_in_cr * cos(angle2),) * 2
    # fifth row
    ham[4][5] = t_in_cr * cos(angle2)
    # tack a minus sign on so the hopping processes lower energy
    return ham


def hamiltonian_cr_cr_3(kx, ky, kz, t_in_cr, t_z_cr, t_out_cr):
    ham = zeros((6, 6), dtype='complex')
    # first row
    ham[0][1] = t_in_cr * fourier(kx, ky, kz, -primitive_vectors[0], angle2)
    ham[0][2] = t_in_cr * fourier(kx, ky, kz, primitive_vectors[1], angle2)
    ham[0][3] = t_z_cr * fourier(kx, ky, kz, primitive_vectors[2], angle1)
    ham[0][4] = t_out_cr * fourier(kx, ky, kz, -primitive_vectors[0] + primitive_vectors[2], angle3)
    ham[0][5] = t_out_cr * fourier(kx, ky, kz, primitive_vectors[1] + primitive_vectors[2], angle3)
    # second row
    ham[1][2] = t_in_cr * fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1], angle2)
    ham[1][3] = t_out_cr * fourier(kx, ky, kz, (primitive_vectors[0] + primitive_vectors[2]), angle3)
    ham[1][4] = t_z_cr * fourier(kx, ky, kz, primitive_vectors[2], angle1)
    ham[1][5] = t_out_cr * fourier(kx, ky, kz, (primitive_vectors[0] + primitive_vectors[1] + primitive_vectors[2]), angle3)
    # third row
    ham[2][3] = t_out_cr * fourier(kx, ky, kz, -primitive_vectors[1] + primitive_vectors[2], angle3)
    ham[2][4] = t_out_cr * fourier(kx, ky, kz, -primitive_vectors[0] - primitive_vectors[1] + primitive_vectors[2], angle3)
    ham[2][5] = t_z_cr * fourier(kx, ky, kz, primitive_vectors[2], angle1)
    # fourth row
    ham[3][4] = t_in_cr * fourier(kx, ky, kz, -primitive_vectors[0], angle2)
    ham[3][5] = t_in_cr * fourier(kx, ky, kz, primitive_vectors[1], angle2)
    # fifth row
    ham[4][5] = t_in_cr * fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1], angle2)
    # tack a minus sign on so the hopping processes lower energy
    return ham


# Constructs the Cr-Cr sector
def hamiltonian_cr_cr(kx, ky, kz, t_a_cr, t_in_cr, t_z_cr, t_out_cr):
    ham1 = full_hamiltonian(hamiltonian_cr_cr_1(kx, ky, kz, t_a_cr))
    ham2 = full_hamiltonian(hamiltonian_cr_cr_2(kx, ky, kz, t_in_cr, t_z_cr, t_out_cr))
    ham3 = full_hamiltonian(hamiltonian_cr_cr_3(kx, ky, kz, t_in_cr, t_z_cr, t_out_cr))
    slice1 = hstack([ham2, zeros((6, 3), dtype='complex'), zeros((6, 6), dtype='complex')])
    slice2 = hstack([zeros((6, 6), dtype='complex'), ham3, zeros((6, 3), dtype='complex')])
    slice3 = hstack([zeros((3, 6), dtype='complex'), zeros((3, 6), dtype='complex'), ham1])
    return vstack([slice1, slice2, slice3])


# The next 3 definitions construct Fe-Fe subsectors
def hamiltonian_fe_fe_1(kx, ky, kz, t_z_fe, t_out_fe):
    ham = zeros((3, 3), dtype='complex')
    # first row
    ham[0][0] = t_z_fe * fourier(kx, ky, kz, primitive_vectors[2], angle5)
    ham[0][1] = t_out_fe * fourier(kx, ky, kz, primitive_vectors[2], angle6)
    ham[0][2] = t_out_fe * fourier(kx, ky, kz, -primitive_vectors[0] - primitive_vectors[1] + primitive_vectors[2], angle6)
    # second row
    ham[1][0] = t_out_fe * fourier(kx, ky, kz, primitive_vectors[2], angle6)
    ham[1][1] = t_z_fe * fourier(kx, ky, kz, primitive_vectors[2], angle5)
    ham[1][2] = t_out_fe * fourier(kx, ky, kz, -primitive_vectors[0] - primitive_vectors[1] + primitive_vectors[2], angle6)
    # third row
    ham[2][0] = t_out_fe * fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1] + primitive_vectors[2], angle6)
    ham[2][1] = t_out_fe * fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1] + primitive_vectors[2], angle6)
    ham[2][2] = t_z_fe * fourier(kx, ky, kz, primitive_vectors[2], angle5)
    return ham


def hamiltonian_fe_fe_2(kx, ky, kz, t_in_fe):
    ham = zeros((3, 3), dtype='complex')
    ham[0][1] = fourier(kx, ky, kz, -primitive_vectors[0], angle4)
    ham[0][2] = fourier(kx, ky, kz, -primitive_vectors[0], angle4)
    ham[1][2] = cos(angle4)
    return t_in_fe * ham


def hamiltonian_fe_fe_3(kx, ky, kz, t_in_fe):
    ham = zeros((3, 3), dtype='complex')
    ham[0][1] = fourier(kx, ky, kz, primitive_vectors[1], angle4)
    ham[0][2] = cos(angle4)
    ham[1][2] = fourier(kx, ky, kz, -primitive_vectors[1], angle4)
    return t_in_fe * ham


# Construct the Fe-Fe sector
def hamiltonian_fe_fe(kx, ky, kz, t_z_fe, t_out_fe, t_in_fe):
    sl1 = hstack([zeros((3, 3), dtype='complex'), hamiltonian_fe_fe_1(kx, ky, kz, t_z_fe, t_out_fe), zeros((3, 9), dtype='complex')])
    sl2 = hstack([adjoint(hamiltonian_fe_fe_1(kx, ky, kz, t_z_fe, t_out_fe)), zeros((3, 12), dtype='complex')])
    sl3 = hstack([zeros((3, 6), dtype='complex'), full_hamiltonian(hamiltonian_fe_fe_2(kx, ky, kz, t_in_fe)), zeros((3, 6), dtype='complex')])
    sl4 = hstack([zeros((3, 9), dtype='complex'), full_hamiltonian(hamiltonian_fe_fe_3(kx, ky, kz, t_in_fe)), zeros((3, 3), dtype='complex')])
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
    return t_cr_fe * ham


def hamiltonian_cr_fe_2(kx, ky, kz, t_cr_fe):
    ham = zeros((3, 3), dtype='complex')
    # first row
    ham[0][0] = fourier(kx, ky, kz, primitive_vectors[0], ang3)
    ham[0][1] = cos(ang2)
    ham[0][2] = cos(ang1)
    # second row
    ham[1][0] = fourier(kx, ky, kz, primitive_vectors[0], ang1)
    ham[1][1] = cos(ang3)
    ham[1][2] = cos(ang2)
    # third row
    ham[2][0] = fourier(kx, ky, kz, primitive_vectors[0], ang2)
    ham[2][1] = cos(ang1)
    ham[2][2] = cos(ang3)
    return t_cr_fe * ham


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
    return t_cr_fe * ham


def hamiltonian_cr_fe_4(kx, ky, kz, t_cr_fe):
    ham = zeros((3, 3), dtype='complex')
    # first row
    ham[0][0] = cos(ang3)
    ham[0][1] = fourier(kx, ky, kz, primitive_vectors[1], ang2)
    ham[0][2] = cos(ang1)
    # second row
    ham[1][0] = fourier(kx, ky, kz, primitive_vectors[0], ang1)
    ham[1][1] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1], ang3)
    ham[1][2] = fourier(kx, ky, kz, primitive_vectors[0], ang2)
    # third row
    ham[2][0] = fourier(kx, ky, kz, -primitive_vectors[1], ang2)
    ham[2][1] = cos(ang1)
    ham[2][2] = fourier(kx, ky, kz, -primitive_vectors[1], ang3)
    return t_cr_fe * ham


def hamiltonian_cr_fe_5(kx, ky, kz, t_cr_fe_p):
    ham = zeros((3, 3), dtype='complex')
    # first row
    ham[0][0] = cos(ang4)
    ham[0][1] = cos(ang4)
    ham[0][2] = fourier(kx, ky, kz, -primitive_vectors[0] - primitive_vectors[1], ang5)
    # second row
    ham[1][0] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1], ang5)
    ham[1][1] = fourier(kx, ky, kz, primitive_vectors[0] + primitive_vectors[1], ang4)
    ham[1][2] = cos(ang4)
    # third row
    ham[2][0] = fourier(kx, ky, kz, primitive_vectors[0], ang4)
    ham[2][1] = fourier(kx, ky, kz, primitive_vectors[0], ang5)
    ham[2][2] = fourier(kx, ky, kz, -primitive_vectors[1], ang4)
    return t_cr_fe_p * ham


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
    return t_cr_fe_p * ham


# Constructs the Cr-Fe sector
def hamiltonian_cr_fe(kx, ky, kz, t_cr_fe, t_cr_fe_p):
    sl1 = hstack([zeros((3, 6), dtype='complex'), hamiltonian_cr_fe_1(kx, ky, kz, t_cr_fe), zeros((3, 6), dtype='complex')])
    sl2 = hstack([zeros((3, 6), dtype='complex'), hamiltonian_cr_fe_2(kx, ky, kz, t_cr_fe), zeros((3, 6), dtype='complex')])
    sl3 = hstack([zeros((3, 9), dtype='complex'), hamiltonian_cr_fe_3(kx, ky, kz, t_cr_fe), zeros((3, 3), dtype='complex')])
    sl4 = hstack([zeros((3, 9), dtype='complex'), hamiltonian_cr_fe_4(kx, ky, kz, t_cr_fe), zeros((3, 3), dtype='complex')])
    sl5 = hstack([hamiltonian_cr_fe_5(kx, ky, kz, t_cr_fe_p), hamiltonian_cr_fe_6(kx, ky, kz, t_cr_fe_p), zeros((3, 9), dtype='complex')])
    return vstack([sl1, sl2, sl3, sl4, sl5])
