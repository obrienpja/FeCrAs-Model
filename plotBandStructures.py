from realisticModel import *

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


# Defines the
def hamiltonian_cr_cr_k_path(t, tz, tperp, tzp, t2p):
    kx, ky, kz = k_path(t)
    return hamiltonian_cr_cr(kx, ky, kz, tz, tperp, tzp, t2p)


def hamiltonian_cr_cr_eigenvalues(t, tz, tperp, tzp, t2p):
    ham = hamiltonian_cr_cr_k_path(t, tz, tperp, tzp, t2p)
    return sorted(map(lambda x: round(x, 5), eigvals(ham)))


def realistic_data(tz, tperp, tzp, t2p):
    return array([hamiltonian_cr_cr_eigenvalues(t, tz, tperp, tzp, t2p) for t in arange(0, 3, 0.01)])


# Plots the Cr-Cr sector over the Gamma-M-K-Gamma-A-H-K-A (Check this)
def realistic_plot(tz, tperp, tzp, t2p):
    data = array([hamiltonian_cr_cr_eigenvalues(t, tz, tperp, tzp, t2p) for t in arange(0, 3, 0.01)])
    plt.gcf().canvas.set_window_title("Realistic Cr-Cr Band Structure Plot")
    plt.title("Realistic Cr-Cr Band Structure Plot")
    plt.xlabel("Kpoints")
    plt.ylabel("Energy")
    plt.axis([0, 300, -6, 2])
    plt.plot(data)
    plt.show()


def read_data(data_file):
    f = open(data_file, 'r')
    data = []
    for line in f:
        data.append(line.split())
    return data


def convert_data(data):
    for i in range(len(data)):
        for j in range(len(data)):
            data[i][j] = complex(data[i][j].replace("*I", "j"))
    return data

dat = array(convert_data(read_data('/home/solidangle/CrCrData.txt')))

for i in range(len(dat)):
    for j in range(len(dat)):
        if dat[i][j] != hamiltonian_cr_cr(.5, .5, .5, 1, 1, 1, 1)[i][j]:
            print i, j, dat[i][j], hamiltonian_cr_cr(.5, .5, .5, 1, 1, 1, 1)[i][j]

print len(dat), len(hamiltonian_cr_cr(.5, .5, .5, 1, 1, 1, 1))
print dat == hamiltonian_cr_cr(.5, .5, .5, 1, 1, 1, 1)

print hamiltonian_cr_cr(.5, .5, .5, 1, 1, 1, 1) == transpose(conj(hamiltonian_cr_cr(.5, .5, .5, 1, 1, 1, 1)))


