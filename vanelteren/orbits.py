from amuse.lab import *
from sandbox.pelupes.kepler3.interface import Kepler
import numpy
from matplotlib import pyplot

def new_orbiters(number_of_orbiters, mass_primary, mass_orbiters_factor = 1e-8, random = None):
    if random is None:
        random = numpy.random
    
    size = (number_of_orbiters,)
    semi_major_axis = (random.random(size) * 0.5 + 0.5) | nbody_system.length
    eccentricity = (random.random(size) * 0.3)
    mass =  (mass_orbiters_factor + (random.random(size) * mass_orbiters_factor)) * mass_primary
    mu = nbody_system.G * (mass + mass_primary)
    
    velocity_perihelion = numpy.sqrt( mu / semi_major_axis  * ((1.0 + eccentricity)/(1.0 - eccentricity)))
    radius_perihelion = semi_major_axis * (1.0 - eccentricity)
    
    result = Particles(number_of_orbiters)
    result.mass = mass
    result.position = radius_perihelion.reshape((number_of_orbiters,1,)) * [[1.0, 0.0, 0.0]]
    result.velocity = velocity_perihelion.reshape((number_of_orbiters,1,)) * [[0.0, 1.0, 0.0]]
    result.semi_major_axis = semi_major_axis
    result.eccentricity = eccentricity
    return result
    
class Orbiters(object):
    dmass_dt = 0.1 | nbody_system.mass / nbody_system.time
    dt = 0.1 | nbody_system.time
    
    # the mass of ther primary cannot get lower than minimum_primary_mass_factor * initial_mass_primary
    minimum_primary_mass_factor = 0.001 
    
    def __init__(self):
        self.all_orbiters_with_error = Particles()
        self.time = None
        
    def setup_model(self, orbiters_at_t0, mass_primary):
        self.initial_mass_primary = mass_primary
        self.code = Kepler()
        self.code.initialize_code()
        self.code.central_particle.add_particle(Particle(
            mass = mass_primary,
            position = [0,0,0] | nbody_system.length,
            velocity = [0,0,0] | nbody_system.speed
        ))
        self.code.orbiters.add_particles(orbiters_at_t0)
        self.previous_state = self.code.orbiters.copy()
            
    def evolve_model(self, tend):
        if self.time is None :
            self.time = 0.0 * tend
        
        if len(self.code.orbiters) == 0:
            raise Exception("no orbiters in model, cannot evolve")
            
        if self.code.central_particle.mass < (self.minimum_primary_mass_factor * self.initial_mass_primary):
            raise Exception("mass of central particle is too low")
            
            
        while self.time < tend:
            try:
                self.code.evolve_model(self.time)
            except Exception as ex:
                self.handle_exception_in_evolve(ex)
                    
            self.previous_state = self.code.orbiters.copy()
            
            self.perform_mass_loss()
            
            self.time += self.dt
            print "time:", self.time, "mass", self.code.central_particle[0].mass
        
    def perform_mass_loss(self):
        self.code.central_particle.mass -= (self.dt * self.dmass_dt)
        if self.code.central_particle.mass < (self.minimum_primary_mass_factor * self.initial_mass_primary):
            raise Exception("mass of central particle is too low")
                
    def get_model(self):
        return self.code.orbiters.copy(), self.all_orbiters_with_error
        
    def stop(self):
        self.code.stop()
    
    def handle_exception_in_evolve(self, ex):
        print ex
        orbiters_with_error = self.code.orbiters[numpy.asarray([len(i) for i in self.code.orbiters.error_message]) > 0]
        print len(orbiters_with_error), "new orbiter(s) with error"
        # need to use orbiters from previous state, as current might be in error
        orbiters_with_error = orbiters_with_error.get_intersecting_subset_in(self.previous_state)
        orbiters_with_error = self.all_orbiters_with_error.add_particles(orbiters_with_error)
        orbiters_with_error.error_time = self.time
        
        self.code.orbiters.remove_particles(orbiters_with_error)
        print len(self.code.orbiters), "orbiter(s) left"
        if len(self.code.orbiters) == 0:
            raise Exception("no orbiters in model, cannot evolve")

def run():
    mass_primary = 1 | nbody_system.mass
    
    orbiters_at_t0 = new_orbiters(10000, mass_primary)
    code = Orbiters()
    code.setup_model(orbiters_at_t0, mass_primary)
    try:
        code.evolve_model(400 | nbody_system.time)
    except Exception as ex:
        print ex
        
    orbits, errored_orbits = code.get_model()
    code.stop()
    
    
    figure = pyplot.figure()
    subplot = figure.add_subplot(1,2,1)
    if len(orbits) > 0:
        subplot.scatter(orbits.semi_major_axis.value_in(nbody_system.length), orbits.eccentricity, c='b')
    subplot.scatter(orbiters_at_t0.semi_major_axis.value_in(nbody_system.length), orbiters_at_t0.eccentricity, c='r')
    
    subplot = figure.add_subplot(1,2,2)
    if len(errored_orbits) > 0:
        subplot.scatter(errored_orbits.semi_major_axis.value_in(nbody_system.length), errored_orbits.eccentricity, c='r')
    
    pyplot.show()
    
if __name__ == '__main__':
    run()
