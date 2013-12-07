"""
   example code for bridging a gravity solver with a hydrodynamics solver
"""
import numpy
from amuse.lab import *
from amuse.couple import bridge
from amuse import datamodel
from amuse.ext.evrard_test import uniform_unit_sphere
from Kepler_orbiters import Orbiters
from amuse.units import quantities

from starplanet import get_kepler_elements, orbital_period, construct_HD106906b

class Mybridge(bridge.Bridge):

    def kick_codes(self,dt):
        #do not kick, sinply return
        #return
        de = quantities.zero
#        print "LP:", self.codes[0].code.central_particle.position
#        print "O:", self.codes[0].code.orbiters.position
        for x in self.codes:
            if hasattr(x,"kick"):
                de += x.kick(dt)

#        print "Mycode=", self.codes[0].__name__

        local_planets = self.codes[1].code.particles.copy()
        local_planets.position -= self.codes[1].code.particles[0].position
        for pi in range(len(local_planets)):
            key_d2 = self.codes[0].code.orbiters.distances_squared(local_planets[pi:pi+1])
            key_coll = numpy.where(key_d2<(1|units.RSun)**2)[0] 
            if len(key_coll)>0:
                print "N collision=", len(key_coll), key_coll, len(local_planets), local_planets[pi], local_planets[pi+1]
                self.codes[0].code.orbiters.remove_particles(self.codes[0].code.orbiters[key_coll])
                print "N left=", len(self.codes[0].code.orbiters)

        self.kick_energy += de

class planet_gravity_for_Heliocentric_orbits(object):
    """
    copy=copycat(base_class,grav_instance, converter)
    derived system, returns copy of grav instance with
    get_gravity_at_point, get_potential_at_point reimplemented in 
    base_class
    """
    def __init__(self,baseclass, sun, planets, converter):
        self.baseclass=baseclass
        self.converter=converter
        self.planets = planets
        self.central_particle = sun
        #print "initialize planet_gravity_for..."
          
    def get_gravity_at_point(self,radius,x,y,z):
        #print "get_gravity_at_point(self,radius,x,y,z)"
        instance=self.baseclass(self.converter)

        # set planets heliocentric
        parts=self.planets.copy()
#        print "pos=", self.central_particle.position
#        print "pos=", parts.position
        parts.position -= self.central_particle.position
        parts.velocity -= self.central_particle.velocity
        instance.particles.add_particles(parts)

        ax,ay,az=instance.get_gravity_at_point(radius,x,y,z)
        
        instance.stop()
        return ax,ay,az

    def get_potential_at_point(self,radius,x,y,z):
        #print "get_get_potential_at_point(self,radius,x,y,z)"
        instance=self.baseclass(self.converter)

        instance.initialize_code()
        parts=self.planets.copy()
        parts.position -= self.central_particle.position
        parts.velocity -= self.central_particle.velocity
        instance.particles.add_particles(parts)

        phi=instance.get_potential_at_point(radius,x,y,z)
        
        instance.stop()
        return phi

def calculate_disk_mass(disk_particles, starpos, Rmin, Rmax):
    
    inner_particles = disk_particles.select(lambda r: (starpos-r).length()<Rmin,["position"])
    outer_particles = disk_particles.select(lambda r: (starpos-r).length()>=Rmax,["position"])

    particles_in_annulus = disk_particles - inner_particles - outer_particles
    
    print "in/out disk:", len(disk_particles), len(inner_particles), len(outer_particles)
    print "Nannulus:", len(particles_in_annulus), particles_in_annulus.mass.sum()

    return particles_in_annulus.mass.sum(), len(particles_in_annulus)


def starplanetwithdisk(peri, apo, t_end, n_steps, Rinner, Router, Mdisk, Ndisk):
    a = 0.5*(peri+apo)
    ecc = (apo-peri)/(apo+peri)
    print "initial orbital parameters:", a, ecc
    bodies = construct_HD106906b(a, ecc)
    bodies.age = zero
    Pplanet = orbital_period(a, bodies.mass.sum())
    sun = bodies[0]

    converter=nbody_system.nbody_to_si(bodies.mass.sum(), a)
    gravity = Huayno(converter) #, redirection="none")
    gravity.particles.add_particles(bodies)
    gravity.parameters.epsilon_squared = (1|units.AU)**2

    channel_from_gravity = gravity.particles.new_channel_to(bodies)
    channel_from_to_gravity = bodies.new_channel_to(gravity.particles)

    dt = t_end/float(n_steps)

    from amuse.ext.protodisk import ProtoPlanetaryDisk
    hydro_converter=nbody_system.nbody_to_si(bodies.mass.sum(), Rinner)
    planetesimals = ProtoPlanetaryDisk(Ndisk, convert_nbody=hydro_converter, 
                                densitypower=1.5, 
                                Rmin=1.0,
                                Rmax=Router/Rinner,
                                q_out=1.0,
                                discfraction=Mdisk/bodies.mass.sum()).result
    print "Mdisk=", Mdisk, planetesimals.mass.sum()
    bodies.move_to_center()
    com = planetesimals.center_of_mass()
    planetesimals.position += bodies[0].position
    planetesimals.velocity += bodies[0].velocity
#    ism.u = 0 | units.m**2 * units.s**-2 
#    ism.h_smooth= 0.01*a

    Pinnerdiskedge = orbital_period(Rinner, bodies[0].mass)
    print "Periods:", Pinnerdiskedge, Pplanet, orbital_period(Router, bodies[0].mass)

    orbiter = Orbiters(converter)
    orbiter.setup_model(planetesimals, sun)
    orbiter.particles = orbiter.orbiters
    orbiter.dmdt = 0 | units.MSun/units.yr
    orbiter.dt = zero

    channel_from_orb_to_framework = orbiter.orbiters.new_channel_to(planetesimals)
    channel_from_framework_to_orb = planetesimals.new_channel_to(orbiter.orbiters)
    channel_from_cp_to_framework = orbiter.central_particle.new_channel_to(planetesimals)
    channel_from_framework_to_cp = planetesimals.new_channel_to(orbiter.central_particle)

    channel_from_orb_to_framework.copy_attributes(["semi_major_axis", "eccentricity"])

    moving_bodies = ParticlesSuperset([bodies, planetesimals])

    filename = "starplanetdisk.hdf5"

    gravorbiter = Mybridge(use_threading=False)
    gravorbiter.add_system(orbiter, (gravity,) )
    gravorbiter.add_system(gravity, () )
    gravorbiter.timestep = min(dt, 0.1*Pplanet)

    Ek0 = gravity.kinetic_energy
    Ep0 = gravity.potential_energy
    E0_tot = Ek0+Ep0

    time = 0 | units.Myr
    while time < t_end:
        time += dt
        bodies.age = time

        write_set_to_file(moving_bodies, filename, 'hdf5')

        gravorbiter.evolve_model(time)
        channel_from_gravity.copy()
        channel_from_orb_to_framework.copy_attributes(["semi_major_axis", "eccentricity"])
#        channel_from_cp_to_framework.copy()

        a, ecc = get_kepler_elements(gravorbiter.model_time, sun, bodies[1], converter) 

        print "Planet:", time, a, ecc
        for pi in planetesimals:
            print "Planetesimal:", time, pi.semi_major_axis, pi.eccentricity

#        Mdisk = calculate_disk_mass(planetesimals, sun, Rinner, Router)
        print "Diskmass=", time, Mdisk
        print "Planetary orbit=", time, a, ecc
        Ek = gravity.kinetic_energy
        Ep = gravity.potential_energy
        E_tot = Ek+Ep
        print "Energies:", E_tot/E0_tot, Ek, Ep


    gravorbiter.stop()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="int", default = 100,
                      help="number of diagnostics time steps [%default]")
    result.add_option("-N", dest="Ndisk", type="int", default = 100,
                      help="number of particles in the disk[%default]")
    result.add_option("-M", unit=units.MJupiter,
                      dest="Mdisk", type="float", default = 3.86984415126e-05|units.MJupiter,
                      help="Mass of the disk [%default]")
    result.add_option("-r", unit=units.AU,
                      dest="Rinner", type="float", default = 50|units.AU,
                      help="inner edge of the disk [%default]")
    result.add_option("-R", unit=units.AU,
                      dest="Router", type="float", default = 120|units.AU,
                      help="outer edge of the disk [%default]")
    result.add_option("--peri", unit=units.AU,
                      dest="peri", type="float", default = 50|units.AU,
                      help="pericenter of planet orbit [%default]")
    result.add_option("--apo", unit=units.AU,
                      dest="apo", type="float", default = 650|units.AU,
                      help="apocenter of planet orbit [%default]")
    result.add_option("-t", unit=units.yr, 
                      dest="t_end", type="float", default = 5000|units.yr,
                      help="end time of the simulation [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.AU, units.yr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")
    o, arguments  = new_option_parser().parse_args()
    starplanetwithdisk(**o.__dict__)

