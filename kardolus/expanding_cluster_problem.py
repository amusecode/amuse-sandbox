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
from amuse.support.codes.core import is_mpd_running

from os import system

from amuse.ext.plummer import MakePlummerModel
from amuse.ext.salpeter import SalpeterIMF

system("echo \" \" > simple_ToverV.dat")

file = open('simple_ToverV.dat', 'a')

class SalpeterIMF(object):
    def __init__(self, mass_min = 0.1 | units.MSun, mass_max = 125 | units.MSun, alpha = -2.35):
        self.mass_min = mass_min.as_quantity_in(units.MSun)
        self.mass_max = mass_max.as_quantity_in(units.MSun)
        self.alpha = alpha
        self.random = random.Random()
        self.random.seed()
    
    def mass_mean(self):
        alpha1 = self.alpha + 1
        alpha2 = self.alpha + 2
        l1 = pow(self.mass_min.value_in(units.MSun), alpha1)
        l2 = pow(self.mass_min.value_in(units.MSun), alpha2)
        u1 = pow(self.mass_max.value_in(units.MSun), alpha1)
        u2 = pow(self.mass_max.value_in(units.MSun), alpha2)
        return ((u2 - l2) * alpha1) / ((u1 - l1) * alpha2) | units.MSun
        
    def mass(self, random_number):
        alpha1 = self.alpha + 1
        factor = (pow(self.mass_max.value_in(units.MSun) / self.mass_min.value_in(units.MSun), alpha1) - 1.0)
        return self.mass_min.value_in(units.MSun) * (pow(1 + (factor * random_number), 1.0 / alpha1)) | units.MSun
        
    def next_mass(self):
        return self.mass(self.random.random())
        
    def next_set(self, number_of_stars):
        set_of_masses = numpy.zeros(number_of_stars)
        total_mass = 0.0 | units.MSun
        for i in range(number_of_stars):
            mass = self.next_mass()
            set_of_masses[i] = mass.value_in(units.MSun)
            total_mass += mass
        return (total_mass, units.MSun.new_quantity(set_of_masses))
        
class SalpeterIMFTests(unittest.TestCase):
    def test1(self):
        instance = SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual(instance.mass_mean().value_in(units.MSun), 0.351, 3)

    def test2(self):
        instance = SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual(instance.mass(1.0).value_in(units.MSun), 100, 3)
        self.assertAlmostEqual(instance.mass(0).value_in(units.MSun), 0.1, 3)
       
    def test3(self):
        instance = SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        n = 10000
        total_mass, set_of_masses = instance.next_set(10000)
        mean = total_mass.value_in(units.MSun) / float(n)
        exact_mean = instance.mass_mean().value_in(units.MSun)
        self.assertTrue(abs(mean - exact_mean) < 0.1)
        
    def test4(self):
        instance = SalpeterIMF(0.1 | units.MSun, 125 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual( 1.0 / instance.mass_mean().value_in(units.MSun), 2.8253, 4)
       

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
    
    particles = MakePlummerModel(number_of_stars, convert_nbody).result;

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
