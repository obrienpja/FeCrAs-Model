from plotBandStructures import *
from random import *

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
font = {'family' : 'normal', 'weight' : 'bold', 'size' : 22}
matplotlib.rc('font', **font)

for i in range(100):
    plt.clf()
    plt.title("OOHM Band Structure for FeCrAs")
    plt.xlabel("Kpoints")
    plt.ylabel("Energy")
    plt.axis([0, 350, -6, 9])
    plt.xticks([0, 50, 100, 150, 200, 250, 300, 350])
    # fig.set_size_inches(15, 10, forward = True)
    non_magnetic_plot(4*random()-2, 4*random()-2, 4*random()-2, 4*random()-2, 4*random()-2, 4*random()-2, 4*random()-2, 4*random()-2, 4*random()-2)
    plt.pause(0.8)
