"""
   Visualization for simple N-body integration.
   Reads particle set from file (nbody.hdf5) and prints subsequent frames.
"""
import sys
import numpy
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser
from time import sleep

def main(filename = "sph.hdf5"):
    pyplot.ion()
    sph_particles = read_set_from_file(filename, 'hdf5')
    lim = max(sph_particles.x).value_in(sph_particles.x.unit)
    for si in sph_particles.history:
        scatter(si.x, si.y)
        xlabel("X")
        ylabel("Y")
        pyplot.xlim(-lim, lim)
        pyplot.ylim(-lim, lim)
#        sleep(0.1)
        pyplot.draw()
        pyplot.cla()
    pyplot.show()

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "sph.hdf5",
                      help="output filename [sph.hdf5]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


