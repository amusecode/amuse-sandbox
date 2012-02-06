"""
   Visualization for simple N-body integration.
   Reads particle set from file (stellar.hdf5) and prints subsequent frames.
"""
import sys
import numpy
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser
from time import sleep

def calculate_effective_temperature(luminosity, radius):
    return ((luminosity/(constants.four_pi_stefan_boltzmann*radius**2))**.25).in_(units.K)

def main(filename = "stellar.hdf5"):
    Tmax = 50000 #| untis.K
    Tmin = 1000 #| untis.K
    Lmax = 1.e+6 #| untis.LSun
    Lmin = 0.1 #| untis.LSun
    pyplot.ion()
    stars = read_set_from_file(filename, 'hdf5')
    stars.add_global_calculated_attribute("temperature", 
                                   calculate_effective_temperature, 
                                   ["luminosity", "radius"])
    m =  10.0*stars.mass/max(stars.mass)
    for si in reversed(list(stars.iter_history())):
        scatter(si.temperature.value_in(units.K), 
                si.luminosity.value_in(units.LSun))#, s=m)
        pyplot.xlabel("T [K]")
        pyplot.ylabel("L [$L_\odot$]")
        pyplot.xlim(Tmax, Tmin)
        pyplot.ylim(Lmin, Lmax)
        pyplot.loglog()
        pyplot.draw()
        pyplot.cla()
    pyplot.show()

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "stellar.hdf5",
                      help="output filename [stellar.hdf5]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


