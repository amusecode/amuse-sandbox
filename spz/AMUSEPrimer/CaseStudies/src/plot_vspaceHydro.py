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

def main(filename = "nbody.hdf5", lim=-1):
    L = lim 
    pyplot.ion()
    storage = store.StoreHDF(filename,"r")
    stars = storage.load()
#    lim = max(stars.x).value_in(stars.x.unit)
    m =  0.05+10.0*stars.mass/max(stars.mass)
    i = 0
    for si in [list(stars.history)[21]]:
        i+=1
        print "i=", i
        si = si.copy()
        si.move_to_center()
        pyplot.figure(figsize=(8,8))
#        pyplot.title("Cluster at t="+str(time))
#        print "time = ", time
        scatter(si[0].vx, si[0].vy, s=100, c='r')
        scatter(si[1].vx, si[1].vy, s=100, c='r')
        scatter(si.vx, si.vy, s=m)
        xlabel("X")
        ylabel("Y")
        if lim>0:
            pyplot.xlim(-lim, lim)
            pyplot.ylim(-lim, lim)
#        lim = 3.e+20
#        pyplot.xlim(-lim, lim)
#        pyplot.ylim(-lim, lim)
        pyplot.draw()
        pyplot.savefig("aaaa.png")
#        sleep(10)
        pyplot.cla()
    pyplot.show()
    sleep(100)

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "nbody.hdf5",
                      help="output filename [nbody.hdf5]")
    result.add_option("-l", dest="lim", type="float", default = -1,
                      help="boxsize")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


