
import numpy
from amuse.ic.plummer import new_plummer_model
from amuse.units import nbody_system
from amuse.community.smalln.interface import SmallN
from amuse.community.kepler.interface import Kepler
from amuse.community.ph4.interface import ph4
from amuse.community.hermite0.interface import Hermite
from amuse import io
from amuse.datamodel import trees
from amuse.datamodel import set_sequential_key_generator

from amuse.units.optparse import OptionParser
from amuse.couple import multiples
from amuse.couple import encounters

import numpy
import os
import logging
import sys

logging.basicConfig(level = logging.DEBUG, stream = sys.stdout)
logging.getLogger("code").setLevel(logging.WARN)
logging.getLogger("amuse.rfi.channel").setLevel(logging.WARN)

class HandleEncounter(encounters.SelectNeighboursByPerturbationMixin, encounters.HandleEncounterWithSmallN):
    
    def __init__(self,
        kepler_code,
        resolve_collision_code,
        interaction_over_code = None,
        G = nbody_system.G
    ):
        encounters.SelectNeighboursByPerturbationMixin.__init__(self)
        encounters.HandleEncounterWithSmallN.__init__(
            self,
            kepler_code,
            resolve_collision_code,
            interaction_over_code,
            G
        )
        
def new_smalln():
    result = SmallN()
    result.parameters.timestep_parameter = 0.1
    return result

def new_kepler():
    code = Kepler()
    code.initialize_code()
    code.parameters.set_defaults()
    return code

def new_ph4():
    code = ph4()
    code.initialize_code()
    code.parameters.set_defaults()
    return code

def new_code_using_ph4(particles):
    gravity_code = new_ph4()
    gravity_code.particles.add_particles(particles)
    gravity_code.stopping_conditions.collision_detection.enable()
    return gravity_code

def new_code_using_multiples(particles):
    gravity_code = new_ph4()
    gravity_code.particles.add_particles(particles)
    gravity_code.stopping_conditions.collision_detection.enable()
    kepler = new_kepler()
    result = multiples.Multiples(gravity_code, new_smalln, kepler)
    result.neighbor_veto = False
    return result
    
def new_code_using_encounters(particles):
    gravity_code = new_ph4()
    encounter_code = HandleEncounter(
        kepler_code =  new_kepler(),
        resolve_collision_code = new_smalln(),
        interaction_over_code = None,
    )
    multiples_code = encounters.Multiples(
        gravity_code = gravity_code,
        handle_encounter_code = encounter_code,
    )
    multiples_code.particles.add_particles(particles)
    multiples_code.commit_particles()
    return multiples_code


def evolve_code(codes, names, dt, endtime,output):
    time = dt * 0
    previous_numbers = [len(codes[0].particles)] * len(names)
    while time < endtime:
        time += dt
        print >> output, "time", time
        reference = None
        new_numbers = []
        for code, name, previous_number_of_particles in zip(codes, names, previous_numbers):
            print "+"*80
            print "+"*10, name
            print "+"*80
            code.evolve_model(time)
            number_of_particles = len(code.particles)
            new_numbers.append(number_of_particles)
            if number_of_particles != previous_number_of_particles:
                print >> output, "for:",name,"multiple found between: ", time - dt, " and: ", time, " number of particles: ", number_of_particles
            if reference is None:
                reference = code.particles.copy()
            else:
                compare = code.particles.copy()
                channel = compare.new_channel_to(reference)
                channel.copy_attribute("x", "x1");
                channel.copy_attribute("y", "y1");
                channel.copy_attribute("z", "z1");
                r2 = (reference.x - reference.x1) ** 2 + (reference.y - reference.y1) ** 2 + (reference.z - reference.z1) ** 2
                
                distances = r2.sqrt()
                i = distances.argmax()
                print >> output, "max key:", reference[i].key
                print >> output, "max distance to reference:", distances.max()
                print >> output, "min distance to reference:", distances.min()
                print >> output, "mean distance to reference:", distances.min()
            print >> output, name, ":number of collisions", code.number_of_collisions
            sys.stdout.flush()
            sys.stderr.flush()
            output.flush()

def parse_arguments():
    parser = OptionParser()
    parser.add_option("-N", dest="number_of_particles", type="int", default=256,
        help="The number of stars in the cluster [%default].")
    parser.add_option("-f", "--radius-factor", dest="radius_factor", type="float", default=1.0,
        help="The factor to change the particle raddi with (by default 1/number_of_particles) [%default].")
    parser.add_option("-d", "--dt", dest="dt", type="float", unit=nbody_system.time, default=0.01 | nbody_system.time,
        help="The delta time to check the number of stars [%default %unit].")
    parser.add_option("-m", dest="number_of_steps", type="int", default = 1000,
        help="The number of steps to take[%default].")
    parser.add_option("-s", dest="seed", type="int", default = 1,
        help="The number of steps to take[%default].")
    parser.add_option("-o", "--output", dest="output_filename", type="string", default = "compare-{number_of_particles}-{seed}-{radius_factor:0>7.2f}.txt",
        help="The name of the output file [%default].")

    options, args = parser.parse_args()
    return options.__dict__

def run(number_of_particles, seed, number_of_steps, dt, radius_factor, output):
    set_sequential_key_generator(1)
    endtime = number_of_steps * dt
    if seed >= 0:
        random = numpy.random.RandomState(seed)
    else:
        random = None
        
    particles = new_plummer_model(number_of_particles, random = random)
    particles.radius = radius_factor * (1.0/number_of_particles) | nbody_system.length

    multiples_code = new_code_using_multiples(particles)
    encounters_code = new_code_using_encounters(particles)
    evolve_code([multiples_code, encounters_code, ],['multiples', 'encounters'], dt, endtime, output)

def print_header(number_of_particles, seed, number_of_steps, dt, radius_factor, output):
    print >> output, "number_of_particles:", number_of_particles
    print >> output, "radius_factor:", radius_factor
    print >> output, "dt:", dt
    print >> output, "number_of_steps:", number_of_steps
    print >> output, "seed:", seed  
    
def main(output, **keyword_arguments):
    print_header(output = output, **keyword_arguments)
    print >> output, "-" * 80
    run(output = output, **keyword_arguments)
    print >> output, "-" * 80
    
def with_output(output_filename, **keyword_arguments):
    output_filename = output_filename.format(**keyword_arguments)
    with open(output_filename, 'w') as stream:
        main(output = stream, **keyword_arguments)
        
if __name__ == "__main__":
    with_output(**parse_arguments())
    
