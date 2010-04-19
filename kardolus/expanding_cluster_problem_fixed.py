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

from amuse.support.units import nbody_system
from amuse.support.units import units

from amuse.legacy.hermite0.interface import Hermite
from amuse.legacy.bhtree.interface import BHTree
from amuse.legacy.phiGRAPE.interface import PhiGRAPEInterface, PhiGRAPE
from amuse.legacy.sse.interface import SSE
from amuse.legacy.support.core import is_mpd_running

from amuse.support.io import store
from os import system

from amuse.ext.plummer import MakePlummerModel
from amuse.ext.salpeter import SalpeterIMF

system("echo \" \" > simple_ToverV.dat")

file = open('simple_ToverV.dat', 'a')
       
def move_particles_to_center_of_mass(particles):
    center_of_mass = particles.center_of_mass()
    center_of_mass_velocity = particles.center_of_mass_velocity()
    
    particles.position = particles.position - center_of_mass
    particles.velocity = particles.velocity - center_of_mass_velocity     


def simulate_small_cluster(number_of_stars, end_time = 40 | units.Myr, name_of_the_figure = "test-2.svg"):
    random.seed()
    
    initial_mass_function = SalpeterIMF()
    total_mass, salpeter_masses = initial_mass_function.next_set(number_of_stars)
    
    convert_nbody = nbody_system.nbody_to_si(total_mass, 1.0 | units.parsec)
    convert_nbody.set_as_default()
    
    print 'end_time:',convert_nbody.to_nbody(end_time)
    print convert_nbody.to_nbody(total_mass)
    
    particles = MakePlummerModel(number_of_stars, convert_nbody).result;

    gravity = PhiGRAPE(convert_nbody,mode="gpu", use_gl="true")
    gravity.parameters.timestep_parameter=0.01
    gravity.parameters.initial_timestep_parameter=0.01
    gravity.parameters.epsilon_squared = 0.0015 | units.parsec ** 2
    
    particles.radius = 0.0 | units.RSun
#    particles.mass = salpeter_masses

    move_particles_to_center_of_mass(particles)

    gravity.particles.add_particles(particles)

    gravity.initialize_particles(0.0)

    gravity.start_viewer()
    
    from_model_to_gravity = particles.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(particles)
   
    time = 0.0 | units.Myr    
#    particles.savepoint(time)

    kineticEnergy = gravity.kinetic_energy.value_in(units.J)
    potentialEnergy = gravity.potential_energy.value_in(units.J)
    ToverV = kineticEnergy / potentialEnergy
    E0=convert_nbody.to_nbody(gravity.kinetic_energy+gravity.potential_energy)
    
    file.write(str(time.value_in(units.Myr)))
    file.write(' ')
    file.write(str(ToverV))
    file.write('\n')

    gravity.parameters.timestep_parameter=0.01
    gravity.parameters.initial_timestep_parameter=0.01

    while time < end_time:
        time += 0.5 * convert_nbody.to_si(1| nbody_system.time)
        gravity.evolve_model(time)
        gravity.synchronize_model()
        from_gravity_to_model.copy()
        print "Evolved model to t    =", str(time)

        kineticEnergy = gravity.kinetic_energy.value_in(units.J)
        potentialEnergy = gravity.potential_energy.value_in(units.J)
        ToverV = kineticEnergy / potentialEnergy
        Et=convert_nbody.to_nbody(gravity.kinetic_energy+gravity.potential_energy)


        print "Kin / Pot             =", ToverV, (Et-E0)/E0 

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
