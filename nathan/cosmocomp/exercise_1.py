import numpy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from amuse.plot import plot, xlabel, ylabel, loglog

from amuse.units import units
from amuse.units import constants


def heating_function(G_0 = 10, eps = 0.05):
    return 10.0**-24 * eps * G_0 | units.erg / units.s

def cooling_function(T, logT = None, a = 3.24, b = 0.170): # x=1e-1
    if logT is None:
        logT = numpy.log10(T.value_in(units.K))
    cooling = 10.0**-21.0 * (10**(-0.1-1.88*(5.23-logT)**4) + 10**(-a-b*abs(4-logT)**3))
    return cooling | units.erg * units.cm**3 / units.s


def plot_cooling_function(cooling_function, figure_name):
    points_per_decade = 20
    logT = numpy.linspace(1.0, 8.0, num = points_per_decade * 7 + 1)
    T = units.K.new_quantity(10**logT)
    
    pyplot.figure(figsize = (12, 10))
    pyplot.gca().set_ylim([1.0e-28, 1.0e-21])
    loglog(T, cooling_function(T))
    xlabel('log(T)')
    ylabel('log($\Lambda$)')
    pyplot.savefig(figure_name)
    print "\nPlot of cooling function was saved to: ", figure_name
    pyplot.close()


if __name__ == '__main__':
    plot_cooling_function(cooling_function, "gerritsen_cooling_function.png")
    
