import math
import numpy
from amuse.lab import *
from optparse import OptionParser

set_printing_strategy("custom",
                      preferred_units = [units.MSun, units.RSun, units.Myr], 
                      precision = 4, prefix = "", separator = " [", suffix = "]")

def calculate_core_mass(star, H_abundance_limit=1.0e-9):
    number_of_zones = star.get_number_of_zones()
    composition = star.get_chemical_abundance_profiles(number_of_zones = number_of_zones)
    # first index where H fraction > 1.0e-9
    index = (composition[0]>H_abundance_limit).nonzero()[0][0] 
    mass = (star.get_cumulative_mass_profile(number_of_zones = number_of_zones) * star.mass)[index]
    return mass[0]

HeWhiteDwarf = 10 | units.stellar_type
def evolve_star_to_core_mass(MZAMS, Mcore, z, H_abundance_limit):
    MZAMS = MZAMS | units.MSun
    Mcore = Mcore | units.MSun

# fix for behaviour where all jobs send their worker to the first machine 
    from socket import gethostname
    #hostname=gethostname(),channel_type="sockets"
    stellar = MESA(hostname=gethostname())
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(Particle(mass=MZAMS))
    stellar.commit_particles()

    Mtcore = 0 | units.MSun
    while Mtcore<Mcore:

        stellar.evolve_model()
        if stellar.particles[0].stellar_type>=HeWhiteDwarf:
            break

        t = stellar.particles[0].age 
        Mt = stellar.particles[0].mass
        Rt = stellar.particles[0].radius

        Mtcore = calculate_core_mass(stellar.particles[0], H_abundance_limit)

    stellar.stop()
    return t, Mt, Rt, Mtcore

def new_option_parser():
    result = OptionParser()
    result.add_option("-H", dest="H_abundance_limit", type="float",
                      default = 1.0e-9,
                      help="Hydrogen abundance limit for core [1.e-9]")
    result.add_option("-M", dest="MZAMS", type="float",default = 3.0,
                      help="Stellar mass [3.0] MSun")
    result.add_option("-m", dest="Mcore", type="float",default = 0.60,
                      help="core mass [0.6] MSun")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="Stellar metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    t, Mt, Rt, Mtcore  = evolve_star_to_core_mass(**o.__dict__)
    print "For ZAMS star of M=", o.MZAMS, "at T=", t, "M=", Mt, "R=", Rt, "Mc=", Mtcore 


