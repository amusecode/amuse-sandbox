"""
   Visualization for simple N-body integration.
   Reads particle set from file (gravity_stellar.hdf5) and prints subsequent frames.
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

def main(filename = "gravity_stellar.hdf5"):
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
        p1 = pyplot.subplot(2,2,1)
        time = si.get_timestamp()
        print time
#        p1.title("Cluster at t="+str(time))
        print "time = ", time
        p1.scatter(si.x.value_in(units.parsec), 
                   si.y.value_in(units.parsec), s=m)
        p1.set_xlabel("X [pc]")
        p1.set_ylabel("Y [pc]")

        p2 = pyplot.subplot(2,2,2)
        p2.scatter(si.temperature.value_in(units.K), 
                si.luminosity.value_in(units.LSun), s=m)
        p2.set_xlabel("T [K]")
        p2.set_ylabel("L [$L_\odot$]")
        p2.set_xlim(Tmax, Tmin)
        p2.set_ylim(Lmin, Lmax)
        p2.loglog()
        pyplot.draw()
        p1.cla()
        p2.cla()
    pyplot.show()

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "gravity_stellar.hdf5",
                      help="output filename [gravity_stellar.hdf5]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


