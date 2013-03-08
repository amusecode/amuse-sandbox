import numpy.random
from amuse.community import *
from amuse.lab import *

#from .interface import eStarsInterface
#from .interface import eStars

from interface import eStarsInterface
from interface import eStars

from matplotlib import pyplot
from amuse.units import units
from amuse.datamodel import Particles
from amuse.ic.brokenimf import new_scalo_mass_distribution
from amuse.ext.particles_with_color import new_particles_with_blackbody_color
from amuse.community.seba.interface import SeBa
from amuse.community.bhtree.interface import BHTree

def new_stellar_evolution(particles):
    stellar_evolution = SeBa()
    stellar_evolution.particles.add_particles(particles)
    return stellar_evolution

def new_gravity(particles, converter):
    gravity = BHTree(converter)
    gravity.particles.add_particles(particles)
    return gravity

if __name__ in ('__main__', '__plot__'):
    number_of_particles = 100
    
    numpy.random.seed(12345)
    masses = new_scalo_mass_distribution(number_of_particles) * 10
    converter = nbody.nbody_to_si(1.0 | units.parsec, masses.sum())
    particles = new_plummer_model(number_of_particles, converter)
    particles.mass = masses
    particles.move_to_center()
    
    gravity = new_gravity(particles, converter)
    stellar_evolution = new_stellar_evolution(particles)
    
    from_gravity_to_local = gravity.particles.new_channel_to(particles)
    from_stellar_evolution_to_local = stellar_evolution.particles.new_channel_to(particles)
    from_stellar_evolution_to_local.copy()
    
    print "adding some color"
    
    particles = new_particles_with_blackbody_color(particles)
    particles.alpha = 1.0
    particles.type = 1
    particles.radius = stellar_evolution.particles.radius.sqrt() * (1e4 | units.parsec).sqrt()
    
    converter = nbody.nbody_to_si(10.0 | units.parsec, masses.sum())
    instance = eStars(converter, channel_type='sockets')#, redirection="none")
    instance.initialize_code()
    instance.particles.add_particles(particles)
    from_local_to_viz = particles.new_channel_to(instance.particles)
    instance.store_view(0|units.Myr)
    
    for i in range(1, 100):
        target_time = i * 0.1 | units.Myr
        gravity.evolve_model(target_time)
        from_gravity_to_local.copy()
        stellar_evolution.evolve_model(target_time) 
        from_stellar_evolution_to_local.copy()
        from_local_to_viz.copy_attributes(["x", "y", "z", "red", "green", "blue"])
        instance.particles.radius = stellar_evolution.particles.radius.sqrt() * (1e4 | units.parsec).sqrt()
        
        instance.store_view(target_time)
    
    instance.stop()
    gravity.stop()
    stellar_evolution.stop()

