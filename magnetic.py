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
    return ham


# Constructs the total Hamiltonian at a given k point
def total_hamiltonian(kx, ky, kz, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe):
    non_mag = kron(identity(2), non_magnetic_hamiltonian(kx, ky, kz, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe))
    return non_mag + magnetic_hamiltonian(mu_val, delta, j_h, sz_cr, sz_fe)


# Constructs the total Hamiltonian at a given BVK point
def total_bvk_hamiltonian(n1, n2, n3, d1, d2, d3, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe):
    return non_magnetic_bvk_hamiltonian(n1, n2, n3, d1, d2, d3, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe) + magnetic_hamiltonian(mu_val, delta, j_h, sz_cr, sz_fe)


# Returns the matrix that diagonalizes the Hamiltonian
def bvk_eigenvectors(n1, n2, n3, d1, d2, d3, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe):
    ham = total_bvk_hamiltonian(n1, n2, n3, d1, d2, d3, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe)
    return eig(ham)[1]


# Plots the realistic model over the Gamma-M-K-Gamma-A-H-L-A (Check this)
def realistic_plot(tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe):
    data = array([total_hamiltonian_eigenvalues(t, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe) for t in arange(0, 3, 0.01)])
    plt.gcf().canvas.set_window_title("Realistic Band Structure Plot")
    plt.title("Realistic Band Structure Plot")
    plt.xlabel("Kpoints")
    plt.ylabel("Energy")
    plt.axis([0, 300, -2, 2])
    plt.plot(data)
    # plt.show()


def muRealistic(mu_Val):
    return total_hamiltonian_eigenvalues(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, mu_Val, .5, 1, 1, 1)


def total_hamiltonian_k_path(t, tz, tperp, tzp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe):
    kx, ky, kz = k_path(t)
    return total_hamiltonian(kx, ky, kz, tz, tperp, tzp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe)


def total_hamiltonian_eigenvalues(t, tz, tperp, tzp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe):
    ham = total_hamiltonian_k_path(t, tz, tperp, tzp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe)
    return sorted(map(lambda x: round(x, 5), eigvals(ham)))


# def non_magnetic_plot(tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, w, h):
#     font = {'family' : 'normal', 'weight' : 'bold', 'size' : 22}
#     matplotlib.rc('font', **font)
#     fig = plt.figure()
#     fig.set_size_inches(w, h, forward = True)
#     data = array([non_magnetic_hamiltonian_eigenvalues(t, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe) for t in arange(0, 6, 0.02)])
#     # plt.gcf().canvas.set_window_title("Realistic Total Band Structure Plot")
#     plt.title("OOHM Band Structure for FeCrAs")
#     plt.xlabel("Kpoints")
#     plt.ylabel("Energy")
#     plt.axis([0, 300, -6, 9])
#     plt.plot(data)
#     plt.show()


# Plots the realistic model over the Gamma-M-K-Gamma-A-H-L-A (Check this)
def realistic_plot(tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe):
    data = array([total_hamiltonian_eigenvalues(t, tz, tperp, txp, t2p, t_cr_fe, t_cr_fe_p, tz_fe, tz_fe_p, tperp_fe, mu_val, delta, j_h, sz_cr, sz_fe) for t in arange(0, 3, 0.01)])
    plt.gcf().canvas.set_window_title("Realistic Cr-Cr Band Structure Plot")
    plt.title("Realistic Cr-Cr Band Structure Plot")
    plt.xlabel("Kpoints")
    plt.ylabel("Energy")
    plt.axis([0, 300, -10, 10])
    plt.plot(data)
    plt.show()


# Convert data from Mathematica's complex form to Python's complex form
def convert_data(data):
    for i in range(len(data)):
        for j in range(len(data)):
            data[i][j] = complex(data[i][j].replace("*I", "j"))
    return data
