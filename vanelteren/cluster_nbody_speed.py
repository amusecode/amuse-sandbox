import sys
import numpy 
import random
import collections
import os
import time

from amuse.support.units import nbody_system
from amuse.support.units import units

from amuse.legacy.hermite0.interface import Hermite
from amuse.legacy.bhtree.interface import BHTree
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from amuse.support.legacy.core import is_mpd_running

from amuse.support.io import store

from amuse.ext.plummer import MakePlummerModel
from amuse.support.data import particle_attributes

def print_log(time, gravity, particles, total_energy_at_t0):
    kinetic_energy = gravity.kinetic_energy
    potential_energy = gravity.potential_energy
    total_energy_at_this_time = kinetic_energy + potential_energy
    print "T     : ", str(time)
    print "ERROR : ", (total_energy_at_this_time - total_energy_at_t0) / total_energy_at_t0
    print "KE    : ", kinetic_energy 
    print "PE    : ", potential_energy
    print "KE/PE : ", kinetic_energy / potential_energy
    print "KE+PE : ", kinetic_energy + potential_energy
    
def simulate_small_cluster(number_of_stars, end_time = 40 | nbody_system.time, number_of_workers = 1):
    particles = MakePlummerModel(number_of_stars).result
    particles.scale_to_standard()
    particles.radius = 0.0 | nbody_system.length
   
    gravity = PhiGRAPE(
        mode = "gpu",
        #number_of_workers = number_of_workers, 
        #debugger = "xterm",
    )
    gravity.initialize_code()
    gravity.parameters.initialize_gpu_once = 1
    
    gravity.particles.add_particles(particles)
    from_gravity_to_model = gravity.particles.new_channel_to(particles)
    
    total_energy_at_t0 = gravity.kinetic_energy + gravity.potential_energy
    
    print "evolving the model until t = " + str(end_time)
    print "gravity evolve step starting"
    for t in end_time.unit.new_quantity(numpy.linspace(0, end_time.number, 1000)):
        gravity.evolve_model(t)
        print t
    print "gravity evolve step done"
        
    from_gravity_to_model.copy()
        
    print_log(time, gravity, particles, total_energy_at_t0)
   
    gravity.stop()
    
    
if __name__ == '__main__':
    result = []
    for number_of_particles in [1024, 2048]: # 256, 512, 1024]:
        for number_of_workers in [1]:#, 2, 3]:
            for ignore in range(1):
                numpy.random.seed(6)
                t0 = time.time()
                simulate_small_cluster(
                    number_of_particles,
                    100.0 | nbody_system.time,
                    number_of_workers
                )
                t1 = time.time()
                measurement = (number_of_particles, number_of_workers, t1-t0)
                print measurement
                result.append(measurement)
        
    print repr(result)
