"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
import numpy

from amuse.lab import *
from optparse import OptionParser
from matplotlib import pyplot

from starplanet import get_kepler_elements 
from starplanetwithdisk import calculate_disk_mass


def plot_starplanet(f1):

    converter=nbody_system.nbody_to_si(1|units.MSun, 100|units.AU)
    pyplot.figure(figsize=(12,12))
    bodies = read_set_from_file(f1, 'hdf5')
    time = zero
    Rinner = 20|units.AU
    Router = 120|units.AU
    t = numpy.array([]) 
    M = numpy.array([]) 
    a = numpy.array([]) 
    r = numpy.array([]) 
    ecc = numpy.array([]) 
    for bi in bodies.history:
        ai, ei = get_kepler_elements(time, bi[0], bi[1], converter) 
        Mdisk, Ndisk = calculate_disk_mass(bi[3:], bi[0].position, Rinner, Router)
        t = numpy.append(t, bi[0].age.value_in(units.yr))
        M = numpy.append(M, Ndisk)
        a = numpy.append(a, ai.value_in(units.AU))
        r = numpy.append(r, bi[1].position.length().value_in(units.AU))
        ecc = numpy.append(ecc, ei)

    Mnorm = float(M[0])
    M = M/Mnorm

    pyplot.subplot(2,2,1)
    pyplot.plot(t, M, c='r', lw=2)
    pyplot.xlabel('t [yr]')
    pyplot.ylabel('M [MSun]')
    pyplot.subplot(2,2,2)
    pyplot.plot(t, a, c='r', lw=2)
    pyplot.plot(t, r, c='b', lw=1)
    pyplot.xlabel('t [yr]')
    pyplot.ylabel('a [AU]')
    pyplot.subplot(2,2,3)
    pyplot.plot(t, ecc, c='r', lw=2)
    pyplot.xlabel('t [yr]')
    pyplot.ylabel('ecc')
    pyplot.show()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="f1",  default = "HDdata.hdf5",
                      help="first file [%default]")

    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.AU, units.yr], 
                          precision = 10, prefix = "", 
                          separator = " [", suffix = "]")

    o, arguments  = new_option_parser().parse_args()
    plot_starplanet(**o.__dict__)
