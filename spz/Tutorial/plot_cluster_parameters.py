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
import operator

MassFraction = [0.005, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] \
    | units.none

def distance_sq(stars) :
    return stars.x**2  + stars.y**2  +stars.z**2

def LagrangianRadii(stars, verbose = 0, massf = MassFraction) :

    com = stars.center_of_mass()
    stars.position = stars.position - com

    vcom = stars.center_of_mass_velocity()
    stars.velocity = stars.velocity - vcom

    n_stars = len(stars)
    # Next line is a potential performance bottleneck, becasue
    # the for loop is not in numpy but explicit
    # old but slow: d2 = numpy.array([distance_sq(star) for star in stars])
    d2 = distance_sq(stars)
    m = stars.mass/stars.mass.sum()
    d2m = zip(d2, m)
    d2m.sort(key = operator.itemgetter(0))

    iL = 0
    mt = 0 | units.none
    Lagrad = []
    for d2i,mi in d2m :
        mt += mi
        while mt >= massf[iL] :
            Lagrad.append(d2i.sqrt())
            # Previous line is preferable above here below, beacuse
            #  the former preserves the units
            #  Lagrad.append(sqrt(d2[ni].value_in(nbody_system.length**2)))
            if verbose==1 :
                print "Lagrangian Radius M= ", mt, \
                      "(iL=", iL, ") at d= ", Lagrad[-1]
            iL += 1
            if iL >= len(massf) :
                break
    return Lagrad

def main(filename = "nbody.hdf5"):
    stars = read_set_from_file(filename, 'hdf5')
    m =  10.0*stars.mass/max(stars.mass)
    time = vector()
    rl01p = AdaptingVectorQuantity() 
    rl25p = AdaptingVectorQuantity() 
    rl50p = AdaptingVectorQuantity() 
    rl75p = AdaptingVectorQuantity() 
    rl90p = AdaptingVectorQuantity() 
    mass = AdaptingVectorQuantity() 
    N = []
    xcom = AdaptingVectorQuantity() 
    ycom = AdaptingVectorQuantity() 
    # make reverse(list(xx)) to reverse the list..
    for si in stars.history:
        time.append(si.get_timestamp())
        si = si.copy()
        com = si.center_of_mass()
        xcom.append(com.x)
        ycom.append(com.y)
        mass.append(si.mass.sum())
        N.append(len(si))
        Rl = LagrangianRadii(si)
        rl01p.append(Rl[1])
        rl25p.append(Rl[5])
        rl50p.append(Rl[6])
        rl75p.append(Rl[7])
        rl90p.append(Rl[8])

    tunit = time.unit
    time = time.value_in(time.unit)
    ax = pyplot.subplot(2,2,1)
    ax.plot(time, mass.value_in(mass.unit))
    ax.set_title("mass")
    ax.set_xlabel("t ["+str(tunit)+"]")
    ax.set_ylabel("M ["+str(mass.unit)+"]")

    ax = pyplot.subplot(2,2,2)
    ax.plot(time, N)
    ax.set_title("N")
    ax.set_xlabel("t ["+str(tunit)+"]")
    ax.set_ylabel("N")

    runit = rl01p.unit
    ax = pyplot.subplot(2,2,3)
    ax.plot(time, rl01p.value_in(runit))
    ax.plot(time, rl25p.value_in(runit))
    ax.plot(time, rl50p.value_in(runit))
    ax.plot(time, rl75p.value_in(runit))
    ax.plot(time, rl90p.value_in(runit))
    ax.set_title("Lagrangian Radii")
    ax.semilogy()
    ax.set_xlabel("t ["+str(tunit)+"]")
    ax.set_ylabel("R ["+str(runit)+"]")

    ax = pyplot.subplot(2,2,4)
    ax.plot(xcom.value_in(runit), ycom.value_in(runit))
    ax.set_title("ceneter of mass")
    ax.set_xlabel("x ["+str(runit)+"]")
    ax.set_ylabel("y ["+str(runit)+"]")

    pyplot.show()

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "nbody.hdf5",
                      help="output filename [nbody.hdf5]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


