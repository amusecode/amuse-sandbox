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

def main(file1 = "nbody.hdf5", file2 = "galaxy.hdf5"):
    pyplot.ion()
    fclu = store.StoreHDF(file1,"r")
    cluster = fclu.load()
    print "cload cluster"
    fgal = store.StoreHDF(file2,"r")
    print "load galaxy"
    galaxy = fgal.load()
    print "loaded"
    m =  10.0*cluster.mass/max(cluster.mass)
    for gi, si in zip(galaxy.history, cluster.history):
#        for i, si in enumerate(cluster.history):
#        gi = galaxy[i].history
        time = si.get_timestamp()
        pyplot.title("Cluster at t="+str(time))
        print "time = ", time
        scatter(gi.x, gi.y,s=0.1)
        scatter(si.x, si.y, c='r',s=m)
#        scatter(si.x, si.y, s=m)
        xlabel("X")
        ylabel("Y")
        lim = 1.e+19
        pyplot.xlim(-lim, lim)
        pyplot.ylim(-lim, lim)
        pyplot.draw()
#        sleep(sleep)
        pyplot.cla()
    pyplot.show()

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="file1", default = "nbody.hdf5",
                      help="output file1 [nbody.hdf5]")
    result.add_option("-F", dest="file2", default = "galaxy.hdf5",
                      help="output filename [nbody.hdf5]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


