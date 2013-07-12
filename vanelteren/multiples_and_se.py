

from amuse.lab import *
from amuse.couple import encounters
from amuse.units import quantities
import numpy
import os
import pdb
import signal
import time
import traceback
import sys

def info(type, value, tb):
    traceback.print_exception(type, value, tb)
    pdb.pm()
    
sys.excepthook = info

def signal_handler(signal, frame):
    print 'You pressed Ctrl+C!'
    pdb.set_trace()
        
def new_smalln(converter):
    result = SmallN(converter)
    result.parameters.timestep_parameter = 0.1
    result.parameters.cm_index = 2001
    return result

def new_kepler(converter):
    kepler = Kepler(converter)
    kepler.initialize_code()
    return kepler

def new_binary(mass1, mass2, semi_major_axis,
               eccentricity = 0, keyoffset = 1):
    total_mass = mass1 + mass2
    mass_fraction_particle_1 = mass1 / (total_mass)
    
#    binary = Particles(keys=range(keyoffset, keyoffset+2))
    binary = Particles(2)
    binary[0].mass = mass1
    binary[1].mass = mass2
    
    mu = constants.G * total_mass
    
    velocity_perihelion = numpy.sqrt( mu / semi_major_axis  * ((1.0 + eccentricity)/(1.0 - eccentricity)))
    radius_perihelion = semi_major_axis * (1.0 - eccentricity)
    print velocity_perihelion
    
    binary[0].position = ((1.0 - mass_fraction_particle_1) * radius_perihelion * [1.0,0.0,0.0])
    binary[1].position = -(mass_fraction_particle_1 * radius_perihelion * [1.0,0.0,0.0])
    
    binary[0].velocity = ((1.0 - mass_fraction_particle_1) * velocity_perihelion * [0.0,1.0,0.0])
    binary[1].velocity = -(mass_fraction_particle_1 * velocity_perihelion * [0.0,1.0,0.0])

    return binary

def kira(tend, N, R):
    PID = os.getpid()
    signal.signal(signal.SIGINT, signal_handler)
    print "You can interrupt this process with:  kill -INT",PID
    time.sleep(5)
    
    set_printing_strategy("custom", 
                      preferred_units = [units.MSun, units.parsec, units.Myr], 
                      precision = 4, prefix = "", 
                      separator = " [", suffix = "]")


#    code = ph4(redirection="none")
    mass = new_salpeter_mass_distribution(N)
    converter = nbody_system.nbody_to_si(mass.sum(), R)
    code = Hermite(converter)
    stars = new_plummer_model(N, convert_nbody=converter)
    stars.mass = mass
    stars.radius = 0.001/len(stars) | R.unit
    stars.velocity *= 0.
    stars[0].mass = 100 | units.MSun
    stars[0].position *= 0
    
    mp = ms = stars[0].mass/2.
    a = 0.5*stars[0].radius
    nb = new_binary(mp, ms, a) 
    print "mas=", nb.mass, nb.key
    binary_particle = stars[0].copy()
    binary_particle.child1 = nb[0]
    binary_particle.child2 = nb[1]
    
    encounter_code = encounters.HandleEncounter(
        kepler_code =  new_kepler(converter),
        resolve_collision_code = new_smalln(converter),
        interaction_over_code = None,
        G=constants.G
        )
    multiples_code = encounters.Multiples(
        gravity_code = code,
        handle_encounter_code = encounter_code,
        G=constants.G
        )
    multiples_code.particles.add_particles(stars)
    multiples_code.singles_in_binaries.add_particles(nb)
    multiples_code.binaries.add_particle(binary_particle)
    multiples_code.commit_particles()
    
    t = quantities.linspace(0*tend, tend, 11)
    for ti in t:
        multiples_code.evolve_model(ti)
        print "at t=", multiples_code.model_time, "N-multiples:", len(multiples_code.multiples)
        # stellar.evolve_mode(ti)
        # binaries.evolve_mode(ti)
        # channel = binaries.particles()...
        # channel.copy
        print "Lagrangian radii:", multiples_code.all_singles.LagrangianRadii(converter)
        print "Lagrangian radii:", multiples_code.particles.LagrangianRadii(converter)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-t", unit=units.Myr,
                      dest="tend",type="float",default=1.|units.Myr)
    result.add_option("-R", unit=units.parsec,
                      dest="R",type="float",default=1|units.parsec)
    result.add_option("-N", 
                      dest="N",type="float",default=100)
    return result

if __name__ == "__main__":
  options, arguments  = new_option_parser().parse_args()
  kira(options.tend, options.N, options.R)
