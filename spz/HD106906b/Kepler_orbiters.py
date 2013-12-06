from amuse.lab import *
from Kepler.interface import Kepler 
from amuse.community.kepler.interface import Kepler as AKepler
import numpy
from matplotlib import pyplot
from amuse.units.quantities import as_vector_quantity

local_kepler_converter=nbody_system.nbody_to_si(1.0|units.MSun, 1|units.AU)
Johannes_Kepler = AKepler(local_kepler_converter, redirection = "none")
Johannes_Kepler.initialize_code()

def orbital_period(a, Mtot) :
    return 2*numpy.pi*(a**3/(constants.G*Mtot)).sqrt()

def get_component_binary_elements(comp1, comp2, kepler):
    mass = comp1.mass + comp2.mass
    pos = comp2.position - comp1.position
    vel = comp2.velocity - comp1.velocity
    kepler.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a,e = kepler.get_elements()
    r = kepler.get_separation()
    E,J = kepler.get_integrals()

    return mass,a,e,r,E

def calculate_orbital_elements(primary, secondary, kepler):
    m,a,e,r,E = get_component_binary_elements(primary, secondary, kepler) 
    m0 = primary.mass
    m1 = secondary.mass
    return a, e, m0, m1

def orbital_parameters_for_the_planets(central_particle, other_bodies, kepler=None, verbose=True):
    if kepler:
        local_kepler = kepler
    else:
        local_kepler = Johannes_Kepler
    a = [] | units.AU
    e = []
    m = [] | units.MSun
    name = []
    for bi in other_bodies:
        ai, ei, M, ms = calculate_orbital_elements(central_particle[0],  bi, local_kepler)
        name.append(bi.name)
        a.append(ai)
        e.append(ei)
        m.append(ms)
    if kepler:
        kepler.stop()
    if verbose:
        for i in range(len(a)):
            print "Planet: ", name[i], a[i], e[i], m[i] 
    other_bodies.semi_major_axis = a
    other_bodies.eccentricity = e

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
    minimum_primary_mass_factor = 0.001 
    
    def __init__(self, converter):
        self.all_orbiters_with_error = Particles()
        self.code = Kepler(converter)
        self.code.initialize_code()
        self.time = None
        
    def setup_model(self, orbiters, primary, time=None):
        self.initial_mass_primary = primary.mass
        if time:
            self.time = time
        else:
            self.time = primary.age
        self.code.central_particle.add_particle(primary)
        self.code.orbiters.add_particles(orbiters)
        self.previous_state = self.code.orbiters.copy()

    @property
    def central_particle(self):
        return self.code.central_particle
    @property
    def orbiters(self):
        return self.code.orbiters

    def evolve_model(self, tend):

        self.orbiters_with_error = Particles()

#        self.dt = min(self.dt, tend-self.time)
        self.dt = tend-self.time
        #print "evolve model:", self.dt, tend, self.time

        if self.time is None :
            self.time = 0.0 * tend
        
        if len(self.code.orbiters) == 0:
            raise Exception("no orbiters in model, cannot evolve")
            
        if self.code.central_particle.mass < (self.minimum_primary_mass_factor * self.initial_mass_primary):
            raise Exception("mass of central particle is too low")
            
            
        while self.time < tend-1.e-5*self.dt:
            try:
                self.code.evolve_model(self.time)
            except Exception as ex:
                self.handle_exception_in_evolve(ex)
                    
            self.previous_state = self.code.orbiters.copy()
            
            self.perform_mass_loss()
            
            self.time += self.dt
#            self.dt = min(self.dt, tend-self.time)
            #print "time:", self.time, "mass", self.code.central_particle[0].mass, self.dt
        
    def perform_mass_loss(self):
        self.code.central_particle.mass += (self.dt * self.dmdt)
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
        
        self.orbiters_with_error = orbiters_with_error
        #self.code.orbiters.remove_particles(orbiters_with_error)
        #print len(self.code.orbiters), "orbiter(s) left"
        """
        if len(orbiters_with_error)>0:
            print "Orbiters with Errors N= ", len(orbiters_with_error)
        self.all_merged_orbiters = Particles()
        if len(orbiters_with_error)>0:
            print orbiters_with_error
            self.all_merged_orbiters.add_particles(orbiters_with_error)
        """
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

def binary_star():
    bodies = Particles(2)
    star = bodies[0]
    planet = bodies[1]

    bodies.age = 0|units.Myr
    star.mass = 1.2|units.MSun
    star.radius = 1|units.RSun
    star.name = "Primary"
    star.position = (0,0,0) | units.AU
    star.velocity = (0,0,0) | units.kms

    planet.mass = 1.0 | units.MSun
    planet.name = "Secondary"
    planet.semi_major_axis = 1|units.AU
    planet.eccentricity = 0.6
    p = planet.semi_major_axis*(1-planet.eccentricity)
    vp =  numpy.sqrt(constants.G*bodies.mass.sum() * (2./p - 1./planet.semi_major_axis))
    planet.position = (1,0,0) * p
    planet.velocity = (0,1,0) * vp

    return bodies

def construct_orbital_orientation(sun, planets, kepler=None):

    if not kepler:
        kepler = Johannes_Kepler

#    R = sun.radius
#    interacting_bodies = planets.select(lambda a, e: a*(1-e)<4*R,["semi_major_axis", "eccentricity"])

    interacting_bodies = planets

    for bi in interacting_bodies:
        kepler.initialize_from_particles(sun+bi)
        lv = kepler.get_longitudinal_unit_vector()
        tv = kepler.get_transverse_unit_vector()
        bi.longitudinal_unit_vector = lv
        bi.transverse_unit_vector = tv #tas_vector_quantity(tv)
#    return interacting_bodies

def reconstruct_posvel(sun, planets, kepler=None):

    if not kepler:
        kepler = Johannes_Kepler

#    R = sun.radius
#    interacting_bodies = planets.select(lambda a, e: a*(1-e)<4*R,["semi_major_axis", "eccentricity"])

    interacting_bodies = planets

    for bi in interacting_bodies:
        mass = sun.mass + bi.mass
        kepler.initialize_from_elements(mass, bi.semi_major_axis, bi.eccentricity)
        lv = bi.longitudinal_unit_vector
        tv = bi.transverse_unit_vector
        kepler.set_longitudinal_unit_vector(lv[0], lv[1], lv[2])
        kepler.set_transverse_unit_vector(tv[0], tv[1], tv[2])
        Porb = orbital_period(bi.semi_major_axis, mass)
        M, A = kepler.get_angles()
        phase = M/(2*numpy.pi) * Porb
        kepler.transform_to_time(phase)

        rpos = as_vector_quantity(kepler.get_separation_vector())
        rvel = as_vector_quantity(kepler.get_velocity_vector())
        bi.position = rpos
        bi.velocity = rvel
    
if __name__ == '__main__':
#    run()

    bodies = binary_star()
    bodies.add_vector_attribute('longitudinal_unit_vector', ('lx','ly','lz')) 
    bodies.add_vector_attribute('transverse_unit_vector', ('tx','ty','tz')) 
    print "0:", bodies
    construct_orbital_orientation(bodies[0:1], bodies[1:])
    print "A:", bodies
    reconstruct_posvel(bodies[0:1], bodies[1:])
    print "B:", bodies
    xxx

    comp = binary_star()
    conv=nbody_system.nbody_to_si(comp.mass.sum(),1|units.AU)
    mass = conv.to_nbody(comp[0].mass + comp[1].mass)
    pos = conv.to_nbody(comp[1].position - comp[0].position)
    vel = conv.to_nbody(comp[1].velocity - comp[0].velocity)
    kepler = AKepler(conv, redirection = "none")
#    kepler = Kepler(conv)
#    kepler = Johannes_Kepler
    kepler.initialize_code()
#    kepler.set_longitudinal_unit_vector(1.0,0.0, 0.0)
#    kepler.set_transverse_unit_vector(0.0, 1.0, 0)
    kepler.initialize_from_particles(comp)
    lv = kepler.get_longitudinal_unit_vector()
    tv = kepler.get_transverse_unit_vector()
    print "lv=", lv
    print "tv=", tv
    xxx

    kepler.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a,e = kepler.get_elements()
    r = kepler.get_separation()
    E,J = kepler.get_integrals()	# per unit reduced mass, note
    lv = kepler.get_longitudinal_unit_vector()
    tv = kepler.get_transverse_unit_vector()
    print "lv=", lv
    print "tv=", tv
    print mass,a,e,r,E

