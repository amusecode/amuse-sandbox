import sys
import unittest
import numpy 
import random
import collections
import os

from amuse.support.units import nbody_system
from amuse.support.units import units

from amuse.legacy.hermite0.interface import Hermite
from amuse.legacy.bhtree.interface import BHTree
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from amuse.legacy.support.core import is_mpd_running

from amuse.support.io import store

from amuse.ext.plummer import MakePlummerModel

def move_particles_to_center_of_mass(particles):
    print  "center of mass:" , particles.center_of_mass()
    print  "center of mass velocity:" , particles.center_of_mass_velocity()
    
    center_of_mass = particles.center_of_mass()
    center_of_mass_velocity = particles.center_of_mass_velocity()
    
    particles.position = particles.position - center_of_mass
    particles.velocity = particles.velocity - center_of_mass_velocity     
    
    print  "center of mass:" , particles.center_of_mass()
    print  "center of mass velocity:" , particles.center_of_mass_velocity()
   
     

def print_log(time, gravity, particles, total_energy_at_t0, total_energy_at_this_time):
    print "T     : " , str(time)
    print "ERROR : " , (total_energy_at_this_time - total_energy_at_t0) / total_energy_at_t0
    print "KE    : " , gravity.kinetic_energy #particles.kinetic_energy(), 
    print "PE    : " , gravity.potential_energy #particles.potential_energy(gravity.parameters.epsilon_squared)
    print "KE/PE : " , gravity.kinetic_energy / gravity.potential_energy
    #print  "center of mass:" , particles.center_of_mass()
    #print  "center of mass velocity:" , particles.center_of_mass_velocity()
    
def simulate_small_cluster(number_of_stars, end_time = 40 | nbody_system.time):
    random.seed()
    
    
    
    particles = MakePlummerModel(number_of_stars, None).result;
   
    gravity = PhiGRAPE(PhiGRAPE.NBODY)
    #gravity = Hermite(Hermite.NBODY)
    
    gravity.setup_module()
    gravity.parameters.epsilon_squared = 0.15 | nbody_system.length ** 2
    gravity.parameters.timestep_parameter = 0.02 | units.none
    gravity.parameters.initial_timestep_parameter = 0.01 | units.none
        
    print "setting masses of the stars"
    particles.radius = 0.0 | nbody_system.length
    print "centering the particles"
    move_particles_to_center_of_mass(particles)
    print  "center of mass:" , particles.center_of_mass()
    
    gravity.particles.add_particles(particles)
    from_model_to_gravity = particles.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(particles)
    
        
    time = 0.0 | end_time.unit
    particles.savepoint(time)
    
    gravity.initialize_particles(0)
    
    total_energy_at_t0 = gravity.kinetic_energy + gravity.potential_energy
    
    print "evolving the model until t = " + str(end_time)
    while time < end_time:
        time += 0.1 | end_time.unit
        
        print "gravity evolve step starting"
        gravity.evolve_model(time)
        print "gravity evolve step done"
        
        #gravity.synchronize_model()
        from_gravity_to_model.copy()

        particles.savepoint(time)  
        
        total_energy_at_this_time = gravity.kinetic_energy + gravity.potential_energy
        print_log(time, gravity, particles, total_energy_at_t0, total_energy_at_this_time)
        
    
    output_file = "nbody.hdf5"
    if os.path.exists(output_file):
        os.remove(output_file)
    storage = store.StoreHDF(output_file)
    storage.store(particles)
   
    del gravity
    
    
        

        
def test_simulate_small_cluster():
    """test_simulate_small_cluster
    This method is found by the testing framework and automatically
    run with all other tests. This method simulates
    a too small cluster, this is done to limit the testing time.
    """
    assert is_mpd_running()
    test_results_path = path_to_test_results.get_path_to_test_results()
    output_file = os.path.join(test_results_path, "test-2.svg")
    simulate_small_cluster(4, 4 | units.Myr, name_of_the_figure = output_file)
    
if __name__ == '__main__':
    simulate_small_cluster(int(sys.argv[1]), int(sys.argv[2]) | nbody_system.time)
