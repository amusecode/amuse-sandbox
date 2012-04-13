import numpy

from matplotlib import pyplot

from amuse.units import nbody_system
from amuse.community.bhtree.interface import BHTree
from amuse.community.gadget2.interface import Gadget2
import logging
import time
from amuse.ic.plummer import new_plummer_model
#logging.basicConfig(level=logging.DEBUG)

smoothing_length = 0.0 | nbody_system.length ** 2


def print_log(time, gravity, particles, total_energy_at_t0):
    kinetic_energy = gravity.kinetic_energy
    potential_energy = gravity.potential_energy
    total_energy_at_this_time = kinetic_energy + potential_energy
    print "time                    : " , time
    print "energy error            : " , (total_energy_at_this_time - total_energy_at_t0) / total_energy_at_t0

    
def simulate_small_cluster(
        number_of_stars = 1000, 
        end_time = 10 | nbody_system.time,
        epsilon_squared = 0.15 | nbody_system.length ** 2
    ):
    numpy.random.seed(1234)
    source = new_plummer_model(number_of_stars)
    
    figure = pyplot.figure()
    subplot = figure.add_subplot(1, 1, 1)
    subplot.scatter(
            source.x.value_in(nbody_system.length),
            source.y.value_in(nbody_system.length),
            s = 1,
            edgecolor = 'blue',
            facecolor = 'blue'
        )
    allparticles = []
    for i, mode in enumerate(("cpu", "g6",)):
        #particles.scale_to_standard()
        print 1
        gravity = BHTree(mode=mode)
        #gravity = Gadget2()
        #gravity.parameters.epsilon_squared = epsilon_squared
        particles = source.copy()
        gravity.particles.add_particles(particles)
        
        from_model_to_gravity = particles.new_channel_to(gravity.particles)
        from_gravity_to_model = gravity.particles.new_channel_to(particles)
        
        total_energy_at_t0 = gravity.kinetic_energy + gravity.potential_energy
        
        positions_at_different_times = []
        positions_at_different_times.append(particles.position)
        #gravity.commit_particles()
        e0 = gravity.kinetic_energy + gravity.potential_energy
        #gravity.evolve_model(1 | nbody_system.time)
        t0 = time.time()
        run_evolve(gravity, end_time)
        t1 = time.time()
        from_gravity_to_model.copy()
        e1 = gravity.kinetic_energy + gravity.potential_energy
        print mode, t1-t0,  (e1 - e0) / e0
        gravity.stop()
        subplot.scatter(
            particles.x.value_in(nbody_system.length),
            particles.y.value_in(nbody_system.length),
            s = 1,
            edgecolor = ['red','green'][i],
            facecolor = ['red','green'][i]
        )
        allparticles.append(particles)
        
    subplot.set_xlim(-4.0,4.0)
    subplot.set_ylim(-4.0,4.0)
    figure.savefig('bhtree.png')
    figure.show()
    difflen = (allparticles[1].position - allparticles[0].position).lengths() / (allparticles[0].position.lengths())
    print "difference in position: %",difflen.mean()*100
    return t1 - t0
    
def run_evolve(gravity, end_time):
    #print "evolving the model form t = {0} to t = {1}".format(gravity.model_time, end_time)
    gravity.evolve_model(end_time)
    

    
    
if __name__ == '__main__':
    for epsilon_squared in (numpy.arange(0.001, 0.002, 0.001) | nbody_system.length ** 2):
        dt = simulate_small_cluster(
            512,
            epsilon_squared = epsilon_squared
        )
        print epsilon_squared / (1|nbody_system.length ** 2), dt
