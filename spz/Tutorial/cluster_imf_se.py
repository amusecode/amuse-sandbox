"""
   Evolve a population of N stars.
   initial mass function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
import sys
import numpy
from amuse.lab import *
from optparse import OptionParser

# import stellar evolution code
from amuse.community.sse.interface import SSE
#from amuse.community.evtwin.interface import EVtwin
    
def main(N=10, t_end=10, dt=1, filename="stellar.hdf5", Mmin=0.1, Mmax= 100, z=0.02):
    t_end = t_end | units.Myr
    dt = dt | t_end.unit
    Mmin = Mmin | units.MSun
    Mmax = Mmax | Mmin.unit

    stellar = SSE()
#    stellar = EVtwin()
    stellar.parameters.metallicity = z

    mZAMS = new_salpeter_mass_distribution(N, Mmin, Mmax)
    bodies = Particles(mass=mZAMS)
    stellar.particles.add_particles(bodies)

    write_set_to_file(stellar.particles, filename, 'hdf5')

    Mtot_init = stellar.particles.mass.sum()
    time = 0.0 | t_end.unit
    while time < t_end:
        time += dt

        stellar.evolve_model(time)
        write_set_to_file(stellar.particles, filename, 'hdf5')

        Mtot = stellar.particles.mass.sum()
        print "T=", time, 
        print "M=", Mtot, "dM[SE]=", Mtot/Mtot_init

    stellar.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-d", dest="dt", type="float", default = 1.0,
                      help="diagnostics time step [1.0] Myr")
    result.add_option("-f", dest="filename", default = "stellar.hdf5",
                      help="output filename [stellar.hdf5]")
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of stars [10]")
    result.add_option("-M", dest="Mmax", type="float",default = 100,
                      help="maximal stellar mass [100] MSun")
    result.add_option("-m", dest="Mmin", type="float",default = 0.1,
                      help="minimal stellar mass [0.1] MSun")
    result.add_option("-t", dest="t_end", type="float", default = 100.0,
                      help="end time of the simulation [100] Myr")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

