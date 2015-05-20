import sys, os
import numpy
import time as pytime
import platform
from numpy import random
from amuse.lab import *
from amuse.couple import bridge
from make_pebble_disk import make_pebble_disk
from amuse.units.optparse import OptionParser
from evolve_disk import *
from make_planets import *
from amuse.units import quantities

data = 10
data * data

MEarth = 5.97219e+24 * units.kg

def new_option_parser():
    result = OptionParser()
    result.add_option("--seed", 
                      dest="seed", type="int", default = 666,
                      help="random number seed [%default]")
    result.add_option("-f", "--file",
                      dest="filename", default = "",
                      help="input filename")
    result.add_option("-t", unit=units.yr,
                      dest="endtime", type="float", default = 2000|units.yr,
                      help="time to run to [%default]")
    result.add_option("--dt", unit=units.yr,
                      dest="dt", type="float", default = 1|units.yr,
                      help="delta time between bridging [%default]")
    result.add_option("-n", "--ndisk",
                      dest="ndisk", type="int", default = 1000,
                      help="Number of pebbels per star [%default]")
    result.add_option("-m", unit=MEarth,
                      dest="disk_mass", type="float", default = 35.0|MEarth,
                      help="disk mass [%default]")
    result.add_option("--inner_disk_radius", unit=units.AU,
                      dest="inner_disk_radius", type="float", default = 15.5 | units.AU,
                      help="inner radius of the pebble disk [%default]")
    result.add_option("--outer_disk_radius", unit=units.AU,
                      dest="outer_disk_radius", type="float", default = 34 | units.AU,
                      help="outer radius of the pebble disk [%default]")
    result.add_option("--code", 
                      dest="code", default = "huayno",
                      help="name of the code to evolve the model with (not used at restart)")
    result.add_option("--new-code", 
                      dest="newcode", default = "",
                      help="name of the code to evolve the model further (used at restart)")
    result.add_option("--particles-kind", 
                      dest="particles_kind", default = "test",
                      help="type of the particles (test, collisionless, gravitational)")
    result.add_option("--no-escapers",
                  action="store_false", dest="with_escapers", default=True,
                  help="don't detect escapers")
    result.add_option("--no-bound-check",
                  action="store_false", dest="check_if_bound_on_escaping", default=True,
                  help="don't check if a particle is bound on escaping")
    result.add_option("--escape_radius", unit=units.AU,
                      dest="escape_radius", type="float", default = 1000000 | units.AU,
                      help="escape radius [%default]")
    return result


class CalculateFieldForParticles(bridge.CalculateFieldForParticles):
    """
    faster algorithm to calculcate the gravity field of a particle set
    assumes:
    1. no epsilon
    2. low number of points, loop in python is removed by doing matrix calculations (needing N x M memory)
    """
    def get_gravity_at_point(self,radius,x,y,z):
        names = ("x","y","z", "mass")
        px, py, pz, mass = self.particles.get_values_in_store(None, names)
        gravity_constant = -self.gravity_constant
        n = len(x)
        newshape =(n, 1)
        x = x.reshape(newshape) 
        y = y.reshape(newshape) 
        z = z.reshape(newshape) 
        dx = x - px
        dy = y - py
        dz = z - pz
        dr_squared = ((dx ** 2) + (dy  ** 2) + (dz ** 2))
        dr_twothird = dr_squared**1.5
        m_div_dr = mass / dr_twothird
        ax = gravity_constant * (m_div_dr*dx).sum(1)
        ay = gravity_constant * (m_div_dr*dy).sum(1)
        az = gravity_constant * (m_div_dr*dz).sum(1)

        
        ax -=  ax[0]
        ay -=  ay[0]
        az -=  az[0]
        #or i,(x, y) in enumerate(zip(ax, ay)):
        #    print i, x.as_quantity_in(units.AU/units.s**2), y.as_quantity_in(units.AU/units.s**2)
        #ddd
        return ax, ay, az
    
    
class GravityCodeInField(bridge.GravityCodeInField):
    """
    faster algorithm to calculcate the gravity for a code
    assumes:
    1. no epsilon
    2. mass of the particles in the code stays constant
    3. no update of the particles in the code (no additions or removals)
    4. use lower level functions: get_values_in_store and set_values_in_store
    """
    def __init__(self, code, field_codes, do_sync=False, verbose=False, radius_is_eps=False, h_smooth_is_eps=False):
        bridge.GravityCodeInField.__init__(self, code, field_codes, do_sync, verbose, radius_is_eps, h_smooth_is_eps)
        self.copy_of_particles = code.particles.copy()
        self.channel_from_code_to_copy = code.particles.new_channel_to(self.copy_of_particles)
        self.channel_from_copy_to_code = self.copy_of_particles.new_channel_to(code.particles)
        self.required_attributes = ['x', 'y', 'z', 'vx', 'vy', 'vz']
        
    def kick_with_field_code(self, particles, field_code, dt):
        
        #positions = particles.position
        names = ('x','y','z', 'vx', 'vy', 'vz')
        x, y, z, vx, vy, vz = particles.get_values_in_store(None, names)
        ax,ay,az=field_code.get_gravity_at_point(
            quantities.zero,
            x,
            y,
            z
        )
        #self.update_velocities(particles, dt, ax, ay, az)
        names = ('vx','vy','vz')
        particles.set_values_in_store(
                None,
                names,
                (vx + dt * ax,
                vy + dt * ay,
                vz + dt * az)
        )
    def update_velocities(self,particles, dt,  ax, ay, az):
        if 1:
            names = ('vx','vy','vz')
            vx, vy, vz = particles.get_values_in_store(None, names)
            particles.set_values_in_store(
                None,
                names,
                (vx + dt * ax,
                vy + dt * ay,
                vz + dt * az)
            )
        else:    
            particles.vx += dt * ax
            particles.vy += dt * ay
            particles.vz += dt * az
        
    def kick(self, dt):
        copy_of_particles = self.copy_of_particles
        self.channel_from_code_to_copy.copy_attributes(self.required_attributes)
        #kinetic_energy_before = copy_of_particles.kinetic_energy()

        for field_code in self.field_codes:
            self.kick_with_field_code(
                copy_of_particles,
                field_code,
                dt
            )

        self.channel_from_copy_to_code.copy_attributes(["vx","vy","vz"])

        #kinetic_energy_after = self.copy_of_particles.kinetic_energy()
        return quantities.zero #kinetic_energy_after - kinetic_energy_before
        
class CreateNiceModel(object):
    
    def __init__(self, number_of_disk_particles = 1000, disk_mass = 35.0 | MEarth, inner_disk_radius = 15.5 | units.AU, outer_disk_radius = 34.4 | units.AU, seed = -1):
        self.number_of_disk_particles = number_of_disk_particles
        self.disk_mass = disk_mass
        self.inner_disk_radius = inner_disk_radius
        self.outer_disk_radius = outer_disk_radius
        self.seed = seed
        if self.seed < 0:
            self.rng = random.RandomState()
        else:
            self.rng = random.RandomState(self.seed)


    def start(self):
        sun = self.create_sun()
        self.result = Particles()
        sun = self.result.add_particle(sun)
        if self.number_of_disk_particles > 0:
            sun.disk_particles = self.create_disk_particles(sun)
        else:
            sun.disk_particles = Particles()
        sun.planets = self.create_planets(sun)
        sun.escaped_disk_particles = Particles()
        sun.escaped_planets = Particles()
        self.result.collection_attributes.inner_disk_radius = self.inner_disk_radius
        self.result.collection_attributes.outer_disk_radius = self.outer_disk_radius
        self.result.collection_attributes.model_time = 0 | units.yr
        return self.result
        
    def create_sun(self):
        sun = Particle()
        sun.mass = 1 | units.MSun
        sun.position = [0,0,0] | units.AU
        sun.velocity = [0,0,0] | units.kms
        sun.radius = 0 | units.AU
        return sun
    
    def create_disk_particles(self, star):
        arange = [self.inner_disk_radius,  self.outer_disk_radius]
        erange = [0., 0.]     
        phi = 0
        theta = 0
        particles = make_pebble_disk(star, self.number_of_disk_particles, arange, erange, phi, theta, self.rng)
        particles.mass = self.disk_mass * 1.0 /len(particles)
        return particles
        

    def create_planets(self, star):
        phi = 0
        theta = 0
        return make_Nice_planets(star, phi, theta, self.rng)
        
        

def potential_energy(pebble, star):
    r = (pebble.position-star.position).length()
    return -constants.G*pebble.mass*star.mass/r

def kinetic_energy(pebble, star):
    return 0.5*pebble.mass * (pebble.velocity-star.velocity).length()**2

def is_bound(pebble, star):
    binding_energy = potential_energy(pebble, star) +  kinetic_energy(pebble, star)
    return binding_energy<zero

class RunNiceModel(object):
    
    def __init__(self, model, dt = None, escape_radius = None, code_name = None, particles_kind = None, with_escapers = None, check_if_bound_on_escaping = None, seed = -1):
        self.model = model
        self.star = self.model[0]
        self.time = self.model.collection_attributes.model_time
        self.time_offset = self.time
        self.channels = []
        self.dt = self._get_value_from_model(dt, self.model, "delta_t", 1 | units.yr)
        self.escape_radius = self._get_value_from_model(escape_radius, self.model, "escape_radius", 1000000 | units.AU)
        self.code_name = self._get_value_from_model(code_name, self.model, "code_name", "hermite")
        self.particles_kind = self._get_value_from_model(particles_kind, self.model, "particles_kind", "collisionless")
        self.with_escapers = self._get_value_from_model(with_escapers, self.model, "with_escapers", True)
        self.check_if_bound_on_escaping = self._get_value_from_model(check_if_bound_on_escaping, self.model, "check_if_bound_on_escaping", True)
        self.seed = seed
        
        

    def _get_value_from_model(self, proposed_value, model, name, default_value):
        if proposed_value is None:
            if hasattr(model.collection_attributes, name):
                return getattr(model.collection_attributes, name)
            else:
                return default_value
        else:
            return proposed_value
            
    def log(self, *args):
        string = ' '.join([str(x) for x in args])
        string += '\n'
        sys.stderr.write(string)
        
    def remove_escaping_particles(self, star, particles, radius, kind = "pebble(s)"):
        escapers = particles[particles.position.lengths()>radius]
        if len(escapers) == 0:
            return escapers
        if self.check_if_bound_on_escaping:
            escapers = escapers[numpy.asarray([is_bound(x, star) for x in escapers])]
        if len(escapers) == 0:
            return escapers
        self.log("removing", len(escapers), " escaping ", kind)
        result = escapers.copy()
        result.escape_time = self.time
        particles.remove_particles(escapers)
        return result
    
    def remove_escaping_pebbles(self, star, radius):
        if len(star.disk_particles) > 0:
            escapers = self.remove_escaping_particles(star, star.disk_particles, radius, "pebble(s)")
            star.escaped_disk_particles.add_particles(escapers)
        
    def remove_escaping_planets(self, star, radius):
        escapers = self.remove_escaping_particles(star, star.planets, radius, "planet(s)")
        star.escaped_planets.add_particles(escapers)
        
    def length_scale(self):
        if len(self.star.disk_particles) == 0:
            return 100 | units.AU
        else:
            return (self.star.disk_particles.position - self.star.position).lengths().max()
        
    def mass_scale(self):
        return self.star.mass
    
    def new_converter(self):
        return nbody_system.nbody_to_si(self.mass_scale(), self.length_scale())
        
    def create_gravity_code(self):
        factory = self.get_code_factory()
        code = factory(self.converter)
        code.particles.add_particles(self.model)
        self.channels.append(code.particles.new_channel_to(self.model))
        
        code.particles.add_particles(self.star.planets)
        self.channels.append(code.particles.new_channel_to(self.star.planets))
        code.particles.move_to_center()
        return code
    
    def get_code_factory(self):
        if self.code_name == "hermite":
            return Hermite
        elif self.code_name == "huayno":
            return Huayno
        elif self.code_name == "mercury":
            return Mercury
            
    
    def create_code(self):
        
        planets_and_star = self.create_gravity_code()
        if len(self.star.disk_particles)  == 0:
            return planets_and_star
        elif self.particles_kind == "gravitational":
            planets_and_star.particles.add_particles(self.star.disk_particles)
            self.channels.append(planets_and_star.particles.new_channel_to(self.star.disk_particles))
            planets_and_star.particles.move_to_center()
            return planets_and_star
        else:
            result = bridge.Bridge(timestep=self.dt, use_threading=True)
            pebble_gravity = advance_without_selfgravity(self.star.disk_particles)
            stellar_gravity = central_point_mass(self.star)
            
            # let the pebbles be kicked by the planets and the sun
            code = GravityCodeInField(pebble_gravity, (planets_and_star,))
            #code = GravityCodeInField(pebble_gravity, (CalculateFieldForParticles(planets_and_star.particles, constants.G),))
            
            result.add_code(code)
            
            if self.particles_kind == "collisionless":
                # let the planets and the sun be kicked by the pebbles
                code = GravityCodeInField(planets_and_star, (CalculateFieldForParticles(self.star.disk_particles, constants.G),))
                result.add_code(code)
            elif self.particles_kind == "test":
                # the pebbles are test particles and do not kick the sun and planets 
                pass
            
            return result
        
    def copy_to_model(self):
        for channel in self.channels:
            channel.copy()
    
    
    def put_parameters_in_model(self):
        self.model.collection_attributes.delta_t = self.dt
        self.model.collection_attributes.code_name = self.code_name
        self.model.collection_attributes.escape_radius = self.escape_radius
        self.model.collection_attributes.particles_kind = self.particles_kind
        self.model.collection_attributes.with_escapers = self.with_escapers
        self.model.collection_attributes.model_time = self.time
        self.model.collection_attributes.check_if_bound_on_escaping = self.check_if_bound_on_escaping
        self.model.collection_attributes.seed = self.seed
        
    

    def get_filename(self):
        return "nice_model-{0}.h5".format(platform.node())
    
    def save(self, append_to_file = True):
        self.put_parameters_in_model()
        write_set_to_file(self.model, self.get_filename(), "amuse", append_to_file=append_to_file, version="2.0")
        
    def run(self, end_time):
        self.converter = self.new_converter()
        self.code = self.create_code()
        istep = 0
        print "starting at time: ", self.time.as_quantity_in(units.yr)
        while self.time < end_time:
            istep += 1
            self.time  += 1000*self.dt
            print "evolving to time: ", self.time.as_quantity_in(units.yr)
            self.code.evolve_model(self.time - self.time_offset)
            for channel in self.channels:
                channel.copy()
            print " evolved to time: ", self.time.in_(units.yr)
            if istep%10 == 0:
                self.save()
            if self.with_escapers and istep%20 == 0:
                self.remove_escaping_pebbles(self.star, self.escape_radius)
                self.remove_escaping_planets(self.star, self.escape_radius)
        
    


def print_model(model):
    print '*' * 50
    print model.collection_attributes
    print '*' * 50
    
if __name__ == "__main__":
    o, arguments  = new_option_parser().parse_args()
    #    random.seed(seed=o.seed)
    
    is_restart = o.filename and os.path.exists(o.filename)
    if is_restart:
        stars = read_set_from_file(o.filename, "amuse", close_file=True).copy()
        if o.newcode:
            stars.collection_attributes.code_name = o.newcode
        
    else:
        uc =  CreateNiceModel(o.ndisk, o.disk_mass, o.inner_disk_radius, o.outer_disk_radius, o.seed)
        uc.start()
        stars = uc.result
        stars.collection_attributes.code_name = o.code
        stars.collection_attributes.particles_kind = o.particles_kind
        stars.collection_attributes.with_escapers = o.with_escapers
        stars.collection_attributes.escape_radius = o.escape_radius
        
        
    
    uc = RunNiceModel(stars, dt = o.dt, seed = o.seed)
    
    
    if not is_restart:
        uc.save(False)
    print_model(uc.model)
    t0 = pytime.time() | units.s
    uc.run(o.endtime)
    t1 = pytime.time() | units.s
    uc.save()
    print_model(uc.model)
    print "total time of this run in seconds: ", t1 - t0
    print "number of days to reach billion years:", \
        (((1.0 | units.Gyr)/uc.time)*(t1 - t0)).as_quantity_in(units.day)

