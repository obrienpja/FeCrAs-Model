from modelFunctions import *
from hoppingSectors import *
from crystalFieldSplitting import *

# Defines origin and special k-points (which ones?)
origin = array([0, 0, 0])
point1 = .5*reciprocal_vectors(*(0,1,2))
point2 = array([.25*norm(reciprocal_vectors(*(0,1,2))), .5*norm(reciprocal_vectors(*(0,1,2))), 0])


# Defines the k path
def k_path(t):
    """kPath definition (had difficulty with the numpy piecewise
    definition of this function. does it matter or is this fine?)"""
    if 0 <= t <= 1:
        return t*point1
    elif 1 < t <= 2:
        return (t - 1) * (point2 - point1) + point1
    elif 2 < t <= 3:
        return (t - 2) * (origin - point2) + point2
    elif 3 < t <= 4:
        return (t - 3) * pi/(2*c) * array([0, 0, 1])
    elif 4 < t <= 5:
        return (t - 4) * point1 + pi/(2*c) * array([0, 0, 1])
    elif 5 < t <= 6:
        return (t - 5) * (point2 - point1) + point1 + pi/(2*c) * array([0, 0, 1])
    elif 6 < t <= 7:
        return (t - 6) * (origin - point2) + point2 + pi/(2*c) * array([0, 0, 1])


# Constructs full non-magnetic Hamiltonian from Cr-Cr, Cr-Fe, Fe-Cr, and Fe-Fe sectors
def non_magnetic_hamiltonian(kx, ky, kz, t_a_cr, t_in_cr, t_z_cr, t_out_cr, t_cr_fe, t_cr_fe_p, t_z_fe, t_out_fe, t_in_fe):
    cfs_cr, cfs_fe = crystal_field_splitting(1000, 1000)
    sl1 = hstack([hamiltonian_cr_cr(kx, ky, kz, t_a_cr, t_in_cr, t_z_cr, t_out_cr) + full_cfs(cfs_cr), hamiltonian_cr_fe(kx, ky, kz, t_cr_fe, t_cr_fe_p)])
    sl2 = hstack([adjoint(hamiltonian_cr_fe(kx, ky, kz, t_cr_fe, t_cr_fe_p)), hamiltonian_fe_fe(kx, ky, kz, t_z_fe, t_out_fe, t_in_fe) + full_cfs(cfs_fe)])
    return vstack([sl1, sl2])


def non_magnetic_hamiltonian_k_path(t, tz, tperp, tzp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe):
    kx, ky, kz = k_path(t)
    return non_magnetic_hamiltonian(kx, ky, kz, tz, tperp, tzp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe)


def non_magnetic_hamiltonian_eigenvalues(t, tz, tperp, tzp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe):
    ham = non_magnetic_hamiltonian_k_path(t, tz, tperp, tzp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe)
    return sorted(map(lambda x: round(x, 5), eigvals(ham)))


def non_magnetic_plot(tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe):
    data = array([non_magnetic_hamiltonian_eigenvalues(t, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe) for t in arange(0, 7, 0.02)])
    # plt.gcf().canvas.set_window_title("Realistic Total Band Structure Plot")
    plt.plot(data)
    # plt.show()
