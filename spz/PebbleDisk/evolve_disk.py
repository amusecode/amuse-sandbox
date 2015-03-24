import sys
import numpy
from numpy import random
from amuse.lab import *
from amuse.couple import bridge
from make_pebble_disk import make_pebble_disk
from amuse.units.optparse import OptionParser

def movie(particles):
    from matplotlib import pyplot
    pyplot.scatter(disk_particles.x.value_in(units.AU), disk_particles.y.value_in(units.AU))
    pyplot.draw()
    pyplot.cla()




        
class central_point_mass(object):
    def __init__(self, particle):
        self.particle = particle

    def get_gravity_at_point(self,eps,x,y,z):
        dx = x-self.particle.x
        dy = y-self.particle.y
        dz = z-self.particle.z

        r2=dx**2 + dy**2 + dz**2
        r=r2**0.5
        fr=constants.G*self.particle.mass/r2
        ax=-fr*dx/r
        ay=-fr*dy/r
        az=-fr*dz/r
        return ax,ay,az

    def get_potential_at_point(self,eps,x,y,z):
        r2=(x-self.particle.x)**2+(y-self.particle.y)**2+(z-self.particle.z)**2
        return constant.G*self.particle.mass/r2

class advance_without_selfgravity(object):
    """
    to advance particles
    """
    def __init__(self, particles, time= 0 |units.yr):
        self._particles = particles
        self.model_time = time

    @property
    def particles(self):
        return self._particles

    def evolve_model(self, t_end):
        dt = t_end - self.model_time
        self._particles.position += self._particles.velocity*dt
        self.model_time= t_end
    
    @property
    def potential_energy(self):
        return quantities.zero
    
    @property 
    def kinetic_energy(self):
        return (0.5*self.particles.mass*self.particles.velocity.lengths()**2).sum()
def new_option_parser():
    result = OptionParser()
    result.add_option("--seed", 
                      dest="seed", type="int", default = 666,
                      help="random number seed [%default]")
    result.add_option("-n", 
                      dest="n", type="int", default = 1000,
                      help="Number of stars [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    random.seed(seed=o.seed)

    arange = [10, 100] | units.AU
    erange = [0., 0.1] 
    Ndisk = o.n
    Mstar = 10|units.MSun
    central_particle = Particle(mass=Mstar, position=(100,0,0)|units.AU, velocity=(0,0,0)|units.kms)

    disk_particles = make_pebble_disk(central_particle, Ndisk, arange, erange)
    disk_particles.mass = 0 | units.MSun

    converter = nbody_system.nbody_to_si(Mstar, arange[-1])

    pebble_gravity = advance_without_selfgravity(disk_particles)
    stellar_gravity = central_point_mass(central_particle)

    dt = 0.1 | units.yr
    gravity = bridge.Bridge(timestep=0.1*dt)
    gravity.add_system(pebble_gravity, (stellar_gravity,), False)
    
    time = zero
    pyplot.ion()
    while time < 100|units.yr:
        time += dt
        gravity.evolve_model(time)
        print time, gravity.model_time
        movie(pebble_gravity.particles)
    pyplot.show()
