import sys
import unittest
import numpy 
import random
import collections
import os
import re

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


from amuse.io import store
from amuse.units import nbody_system
from amuse.units import units
from amuse.community.hermite0.interface import Hermite
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPEInterface, PhiGRAPE
from amuse.community.sse.interface import SSE
from os import system

from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_sphere
from amuse.ic.salpeter import new_salpeter_mass_distribution
system("echo \" \" > simple_ToverV.dat")

file = open('simple_ToverV.dat', 'a')


def move_particles_to_center_of_mass(particles):
    center_of_mass = particles.center_of_mass()
    center_of_mass_velocity = particles.center_of_mass_velocity()
    
    particles.position = particles.position - center_of_mass
    particles.velocity = particles.velocity - center_of_mass_velocity     


def simulate_small_cluster(number_of_stars, end_time = 40 | units.Myr, name_of_the_figure = "test-2.svg"):
    random.seed()
    
    salpeter_masses = new_salpeter_mass_distribution(number_of_stars)
    total_mass = salpeter_masses.sum()
    
    convert_nbody = nbody_system.nbody_to_si(total_mass, 1.0 | units.parsec)
    convert_nbody.set_as_default()
    
    particles = new_plummer_sphere(number_of_stars, convert_nbody);

    gravity = PhiGRAPE(convert_nbody,mode="gpu")
    gravity.initialize_code()
    gravity.parameters.timestep_parameter=0.01
    gravity.parameters.initial_timestep_parameter=0.01 
    gravity.parameters.epsilon_squared = 0.000001 | units.parsec ** 2
     
    stellar_evolution = SSE()
    stellar_evolution.initialize_module_with_default_parameters() 
    
    #print "setting masses of the stars"
    particles.radius = 0.0 | units.RSun
    
    #comment out to get plummer masses
    particles.mass = salpeter_masses
    
    gravity.particles.add_particles(particles)
    gravity.initialize_particles(0.0)

    from_model_to_gravity = particles.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(particles)
   
    time = 0.0 | units.Myr    
    particles.savepoint(time)
 
    move_particles_to_center_of_mass(particles)

    kineticEnergy = gravity.kinetic_energy.value_in(units.J)
    potentialEnergy = gravity.potential_energy.value_in(units.J)
    ToverV = kineticEnergy / potentialEnergy

    file.write(str(time.value_in(units.Myr)))
    file.write(' ')
    file.write(str(ToverV))
    file.write('\n')

    while time < end_time:
        time += 0.25 | units.Myr
        gravity.evolve_model(time)
        from_gravity_to_model.copy()
        print "Evolved model to t    =", str(time)
        print "Evolved model to t    =", str(convert_nbody.to_nbody( time.value_in(units.Myr)| units.Myr))

        kineticEnergy = gravity.kinetic_energy.value_in(units.J)
        potentialEnergy = gravity.potential_energy.value_in(units.J)
        ToverV = kineticEnergy / potentialEnergy

        print "Kin / Pot             =", ToverV
        #print "Particle Mass         =", particles[1].mass

        file.write(str(time.value_in(units.Myr)))
        file.write(' ')
        file.write(str(ToverV))
        file.write('\n')
        
    file.close()      
    
    if os.path.exists('small.hdf5'):
        os.remove('small.hdf5')
    storage = store.StoreHDF("small.hdf5")
    storage.store(particles)
   
    del gravity
    del stellar_evolution
    
        
def test_simulate_small_cluster():
    """test_simulate_small_cluster
    This method is found by the testing framework and automatically
    run with all other tests. This method simulates
    a too small cluster, this is done to limit the testing time.
    """
    assert is_mpd_running()
    simulate_small_cluster(4, 4 | units.Myr)
    
if __name__ == '__main__':
    simulate_small_cluster(int(sys.argv[1]), int(sys.argv[2]) | units.Myr)#, sys.argv[3])
