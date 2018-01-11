from plotBandStructures import *
from random import *

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
font = {'family' : 'normal', 'weight' : 'bold', 'size' : 22}
matplotlib.rc('font', **font)

def ani(ran, num):
    for i in range(num):
        plt.clf()
        plt.title("OOHM Band Structure for FeCrAs")
        plt.xlabel("Kpoints")
        plt.ylabel("Energy")
        plt.axis([0, 350, -5, 8])
        plt.xticks([0, 50, 100, 150, 200, 250, 300, 350],['G','M','K','G','A','L','H','A'])
        # fig.set_size_inches(15, 10, forward = True)
        # non_magnetic_plot(ran*random()-ran/2, ran*random()-ran/2, ran*random()-ran/2, ran*random()-ran/2, ran*random()-ran/2, ran*random()-ran/2, ran*random()-ran/2, ran*random()-ran/2, ran*random()-ran/2)
        non_magnetic_plot(i*10, 1, 1, 1, 1, 1, 1, 1, 1)
        plt.pause(0.8)
