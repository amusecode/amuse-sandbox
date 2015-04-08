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


def remove_escapers(stars, radius):
    for star in stars:
        pebbles = star.disk_particles
        escapers = pebbles[pebbles.position.lengths()>radius]
        if len(escapers) == 0:
            return

        if 0:
            # we no longer check if an escaper is still bound to a star (as checking
            # if a pebble is bound is time consuming
            escapers = escapers[numpy.asarray([x is None for x in escapers.star])]
            if len(escapers)>0:
                pebbles.remove_particles(escapers)
        else:
            print "removing", len(escapers), "escaping pebble(s)"
            pebbles.remove_particles(escapers)
            
def remove_escaping_planets(stars, radius):
    for star in stars:
        planets = star.planets
        escapers = planets[planets.position.lengths()>radius]
        if len(escapers) == 0:
            return
        
        if 0:
            escapers = escapers[numpy.asarray([x is None for x in escapers.star])]
            if len(escapers)>0:
                planets.remove_particles(escapers)
        else:
            print "removing", len(escapers), "escaping planet(s)"
            planets.remove_particles(escapers)

def new_option_parser():
    result = OptionParser()
    result.add_option("--seed", 
                      dest="seed", type="int", default = 666,
                      help="random number seed [%default]")
    result.add_option("-N", 
                      dest="Nstars", type="int", default = 1,
                      help="Number of stars [%default]")
    result.add_option("-f", 
                      dest="filename", default = "",
                      help="input filename")
    result.add_option("-R", unit=units.parsec,
                      dest="Rcluster", type="int", default = 1000|units.AU,
                      help="Cluster size [%default]")
    result.add_option("-t", unit=units.yr,
                      dest="endtime", type="float", default = 2000|units.yr,
                      help="time to run to [%default]")
    result.add_option("-n", 
                      dest="n", type="int", default = 1000,
                      help="Number of pebbels per star [%default]")
    return result


class CalculateFieldForParticles(bridge.CalculateFieldForParticles):
    """
    faster algorithm to calculcate the gravity field of a particle set
    assumes:
    1. no epsilon
    2. low number of points, loop in python is removed by doing matrix calculations (needing N x M memory)
    """
    def get_gravity_at_point(self,radius,x,y,z):
        positions = self.particles.position
        mass = self.particles.mass
        gravity_constant = -self.gravity_constant
        n = len(x)
        newshape =(n, 1)
        x = x.reshape(newshape) 
        y = y.reshape(newshape) 
        z = z.reshape(newshape) 
        dx = x - positions.x
        dy = y - positions.y
        dz = z - positions.z
        dr_squared = ((dx ** 2) + (dy  ** 2) + (dz ** 2))
        dr_twothird = dr_squared**1.5
        ax = gravity_constant * (mass*dx/dr_twothird).sum(1)
        ay = gravity_constant * (mass*dy/dr_twothird).sum(1)
        az = gravity_constant * (mass*dz/dr_twothird).sum(1)

        
        ax -=  ax[0]
        ay -=  ay[0]
        az -=   az[0]
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
        names = ('x','y','z')
        x, y, z = particles.get_values_in_store(None, names)
        ax,ay,az=field_code.get_gravity_at_point(
            quantities.zero,
            x,
            y,
            z
        )
        self.update_velocities(particles, dt, ax, ay, az)

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
        kinetic_energy_before = copy_of_particles.kinetic_energy()

        for field_code in self.field_codes:
            self.kick_with_field_code(
                copy_of_particles,
                field_code,
                dt
            )

        self.channel_from_copy_to_code.copy_attributes(["vx","vy","vz"])

        kinetic_energy_after = self.copy_of_particles.kinetic_energy()
        return kinetic_energy_after - kinetic_energy_before
    
if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    #    random.seed(seed=o.seed)
    
    if o.filename and os.path.exists(o.filename):
        stars = read_set_from_file(o.filename, "amuse", close_file=True).copy()
        Rcl = o.Rcluster
        converter = nbody_system.nbody_to_si(stars.mass.sum(), Rcl)
        
        time = stars.collection_attributes.model_time
        begin_time = time
    else:
        Ncl = o.Nstars
        masses = 1.0 | units.MSun
        Mcl = masses.sum()
        Rcl = o.Rcluster
        converter = nbody_system.nbody_to_si(Mcl, Rcl)
        stars = Particles(1)
        stars.mass = masses
        stars.position = [0,0,0] | units.AU
        stars.velocity = [0,0,0] | units.kms
        stars.radius = 0 | units.AU
        arange = [15.5, 34] | units.AU
        erange = [0., 0.] 
        Ndisk = o.n
        for star in stars:
#            phi = numpy.radians(random.uniform(0, 90, 1)[0])#rotate under x
#            theta = numpy.radians(random.uniform(0, 180, 1)[0]) #rotate under y
            phi = 0
            theta = 0
            star.disk_particles = make_pebble_disk(star, Ndisk, arange, erange, phi, theta)
            star.planets = make_Nice_planets(star, phi, theta)
            MEarth = 5.97219e+24 * units.kg
            Mdisk = 35.0 | MEarth
            star.disk_particles.mass = Mdisk/len(star.disk_particles)
            star.bound_pebbles = Particles()
        time = 0 | units.yr
        begin_time = time

    
#    cluster = Hermite(converter)
#    cluster.parameters.dt_param = 0.0003
#    cluster = Huayno(converter)
#    cluster.parameters.inttype_parameter = 21
    cluster = Hermite(converter)
#    cluster.parameters.begin_time = time
#    from amuse.community.sakura.interface import Sakura
#    cluster = Sakura(converter)
    cluster.particles.add_particles(stars)
    channel_from_cluster = cluster.particles.new_channel_to(stars)
    tcr = Rcl/numpy.sqrt(constants.G*stars.mass.sum()/Rcl)
    print "Tcr=", tcr.in_(units.yr)
    dt = 1|units.yr
    cluster_with_pebble_disks = bridge.Bridge(timestep=1.0 * dt, use_threading=True)

    channels_from_planets = []
    for star in stars:
        cluster.particles.add_particles(star.planets)
        channels_from_planets.append(
            cluster.particles.new_channel_to(star.planets))
    cluster.particles.move_to_center()
    
    for star in stars:
        pebble_gravity = advance_without_selfgravity(star.disk_particles)
        stellar_gravity = central_point_mass(star)
        # let the pebbles be kicked by the planets and the sun
        code = GravityCodeInField(pebble_gravity, (cluster,))
        cluster_with_pebble_disks.add_code(code)
        # let the planets and the sun be kicked by the pebbles
        code = GravityCodeInField(cluster, (CalculateFieldForParticles(star.disk_particles, constants.G),))
        cluster_with_pebble_disks.add_code(code)
        
        # let the planets and the sun NOT be kicked by the pebbles
        #code = GravityCodeInField(cluster, [])
        #cluster_with_pebble_disks.add_code(code)
        
#    cluster_with_pebble_disks.add_system(cluster, ())
    t0 = time
    filename = "nice_model-{0}-{1}.h5".format(platform.node(), int(t0.value_in(1000 * units.yr)))
    write_set_to_file(stars, filename, "amuse", append_to_file=False, version="2.0")
    
    istep = 0
    t0 = pytime.time() | units.s
    print "starting at time: ", time.as_quantity_in(units.yr)
    while time < o.endtime:
        istep += 1
        time += 1000*dt
        print "evolving to time: ", time.as_quantity_in(units.yr)
        cluster_with_pebble_disks.evolve_model(time - begin_time)
        channel_from_cluster.copy()
        for channel in channels_from_planets:
            channel.copy()
        print "evolved to time: ", time.in_(units.yr)
        if istep%10 == 0:
            write_set_to_file(stars, filename, "amuse", append_to_file=True, version="2.0", extra_attributes = {'model_time': time})
        if False and istep%20 == 0:
            remove_escapers(stars, 1000000|units.AU)
            remove_escaping_planets(stars, 1000000|units.AU)
            
    
    t1 = pytime.time() | units.s
    print "total time of this run in seconds: ", t1 - t0
    print "number of days to reach billion years:", \
        (((1.0 | units.Gyr)/time)*(t1 - t0)).as_quantity_in(units.day)
