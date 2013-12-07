"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
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
    Rinner = 50|units.AU
    Router = 120|units.AU
    t = []
    M = [] 
    i = 0
    for bi in bodies.history:
        a, ecc = get_kepler_elements(time, bi[0], bi[1], converter) 
        Mdisk = calculate_disk_mass(bi[3:], bi[0].position, Rinner, Router)
        i+=1
        t.append(i)
        M.append(Mdisk.value_in(units.MJupiter))
        

    pyplot.plot(t, M, c='r', lw=2)
    pyplot.xlabel('t [i]')
    pyplot.ylabel('M [MSun]')
    pyplot.show()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="f1",  default = "HD106906_data.hdf5",
                      help="first file [%default]")

    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.AU, units.yr], 
                          precision = 10, prefix = "", 
                          separator = " [", suffix = "]")

    o, arguments  = new_option_parser().parse_args()
    plot_starplanet(**o.__dict__)
