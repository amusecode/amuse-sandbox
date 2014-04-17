import time
import numpy

from amuse.units import nbody_system
from amuse.units import quantities
from amuse.community.hermite0.interface import Hermite
from amuse.community.kepler.interface import Kepler
from amuse.community.smalln.interface import SmallN
from amuse.couple import encounters
from amuse.ic import plummer

import logging

from amuse.units.optparse import OptionParser

class CountingHandleEncounter(encounters.HandleEncounter):
    
    def __init__(self,
        kepler_code,
        resolve_collision_code,
        interaction_over_code = None,
        G = nbody_system.G
    ):
        encounters.HandleEncounter.__init__(
            self,
            kepler_code,
            resolve_collision_code,
            interaction_over_code,
            G
        )
        self.number_of_encounters = 0

    def execute(self):
        encounters.HandleEncounter.execute(self)
        self.number_of_encounters += 1
        
    def get_number_of_encounters(self):
        return self.number_of_encounters

def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-n", "--number_of_particles", 
        default = 100,
        dest="number_of_particles",
        help="number of particles in the cluster",
        type="int"
    )
    result.add_option(
        "-s", "--seed", 
        default = -1,
        dest="seed",
        help="random seed to use",
        type="int"
    )
    result.add_option(
        "-t", "--time", 
        default = 10 | nbody_system.time,
        dest="end_time",
        help="nbody time to run to",
        unit=nbody_system.time,
        type="float"
    )
    result.add_option(
        "-r", "--radius", 
        default = 0.01 | nbody_system.length,
        dest="radius",
        help="radius to use",
        unit=nbody_system.length,
        type="float"
    )
    return result
    
def run_hermite(number_of_particles, seed, end_time_nbody, radius):
    if seed > 0:
        numpy.random.seed(seed)
    
    scale =  1 | nbody_system.length
    
    model = plummer.new_plummer_model(number_of_particles)
    
    kepler = Kepler()
    kepler.initialize_code()
    
    smalln = SmallN()
    smalln.parameters.timestep_parameter = 0.1
    smalln.parameters.cm_index = 2001
    
        
    encounter_code = CountingHandleEncounter(
        kepler_code =  kepler,
        resolve_collision_code = smalln,
        interaction_over_code = None
    )
    
    gravity_code = Hermite()
    code = encounters.Multiples(
        gravity_code = gravity_code,
        handle_encounter_code = encounter_code
    )
    model.radius = radius
    code.particles.add_particles(model)
    end_time = end_time_nbody
    code.commit_particles()
    e0 = code.get_total_energy()
    dt = end_time / 100.0
    t = dt
    print "starting..."
    try:
        while t < end_time:
            code.evolve_model(t)
            print t, code.get_total_energy(), (code.get_total_energy() - e0) / e0, len(code.multiples)
            particles = code.all_singles
            if 0:
                coreposition,coreradius,coredens=particles.densitycentre_coreradius_coredens()
                print "--", code.model_time.value_in(nbody_system.time),  coreposition.length().value_in(nbody_system.length), coreradius.value_in(nbody_system.length), coredens.value_in(nbody_system.mass / nbody_system.length**3)
            
            radii = particles.LagrangianRadii(mf=[0.1,0.5])
            print "--", code.model_time.value_in(nbody_system.time), radii[0][0].value_in(nbody_system.length), radii[0][1].value_in(nbody_system.length),(code.get_total_energy() - e0) / e0, len(code.multiples), encounter_code.get_number_of_encounters()
            t += dt
    finally:
        print "stopping..."
        gravity_code.stop()
        kepler.stop()
        smalln.stop()
        print "done"
       
if  __name__ == '__main__':
    
    options, arguments = new_option_parser().parse_args()
    
    logging.basicConfig(level=logging.ERROR)
    #encounters.LOG_ENERGY.setLevel(logging.DEBUG)
    run_hermite(options.number_of_particles, options.seed, options.end_time, options.radius)
