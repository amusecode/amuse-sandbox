"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""

from  time import clock
import numpy 
from amuse.lab import *
from optparse import OptionParser
from starplanet import construct_HD106906b
from Kepler.interface import Kepler
from Kepler_orbiters import orbital_period, orbital_parameters_for_the_planets, construct_orbital_orientation, reconstruct_posvel


from amuse.couple import bridge
from amuse.units import quantities

def print_star(star):
    print "Star: t=",  star.age, 
    print "m=", star.mass, 
    print "r=", star.radius 

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

        local_planets = self.codes[1].code.particles.copy()
        local_planets.position -= self.codes[1].code.particles[0].position
        for pi in range(len(local_planets)):
            key_d2 = self.codes[0].code.orbiters.distances_squared(local_planets[pi:pi+1])
            key_coll = numpy.where(key_d2<(1|units.RSun)**2)[0] 
#            print "N=", key_d2, len(key_d2), key_coll, pi, pi+1
            if len(key_coll)>0:
                print "N collision=", len(key_coll), key_coll, local_planets[pi], local_planets[pi+1]
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

def remove_some_particles(remove_bodies, remove_from_sets):

    if len(remove_bodies)>0:
        print "remove particles N=", len(remove_bodies)
        print remove_bodies.key
        for x in remove_from_sets:
            y = x.difference(remove_bodies)
            if len(y)<len(x):
                print "remove particle from set N=", len(x), len(y)
                print x.key
                x.remove_particles(remove_bodies)


def merge_two_particles(gravity, p1, p2, bodies):
    particles_in_encounter = Particles(particles=[p1, p2])
#    particles_in_encounter = particles_in_encounter.get_intersecting_subset_in(bodies)

    com_pos = particles_in_encounter.center_of_mass()
    com_vel = particles_in_encounter.center_of_mass_velocity()
    gravity.particles[0].mass += p2.mass
    gravity.particles[0].position = com_pos
    gravity.particles[0].velocity = com_vel
    gravity.particles.remove_particle(p2)
    gravity.particles.synchronize_to(bodies)

def integrate_planetary_system(sun_and_planets, planetesimals, t_end, n_steps, verbose):

    converter=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
    orbital_parameters_for_the_planets(sun_and_planets[0:1], sun_and_planets[1:], verbose=False)
    orbital_parameters_for_the_planets(sun_and_planets[0:1], planetesimals, verbose=False)

    pd_gravity = Kepler(converter)
    ss_gravity = Huayno(converter)
    bodies = ParticlesSuperset([sun_and_planets, planetesimals])
    bodies.age = 0|units.yr

    ss_gravity.particles.add_particles(sun_and_planets)
    ss_gravity.commit_particles()

    sun = sun_and_planets[0:1]
    pd_gravity.central_particle.add_particles(sun)
    pd_gravity.central_particle.position = sun_and_planets[0].position
    pd_gravity.central_particle.velocity = sun_and_planets[0].velocity
    pd_gravity.orbiters.add_particles(planetesimals)
    pd_gravity.particles = pd_gravity.orbiters
    pd_gravity.commit_particles()

#    print "a, e=", pd_gravity.orbiters.get_eccentricity()

    trun_start = clock() | units.s
    one_step_for_planetary_system(sun_and_planets, planetesimals, 
                                  ss_gravity, pd_gravity, 
                                  t_end, n_steps, verbose)
    trun_stop = clock() | units.s
    print "Performace data: N=", len(sun_and_planets), "n=", len(planetesimals), "t_end=", t_end, "dt=", trun_stop - trun_start

    orbital_parameters_for_the_planets(bodies[0:1], bodies[1:], verbose=False)

    ss_gravity.stop()
    pd_gravity.stop()

def one_step_for_planetary_system(sun_and_planets, planetesimals, 
                                  planets_gravity, planetesimal_gravity, 
                                  t_end, n_steps, verbose):

    bodies = ParticlesSuperset([sun_and_planets, planetesimals])
    bodies.age = 0|units.yr

    channel_from_ssd_to_framework = planets_gravity.particles.new_channel_to(sun_and_planets)
    channel_from_framework_to_ssd = sun_and_planets.new_channel_to(planets_gravity.particles)
    channel_from_pdd_to_framework = planetesimal_gravity.particles.new_channel_to(planetesimals)
    channel_from_framework_to_pdd = planetesimals.new_channel_to(planetesimal_gravity.particles)

    converter=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)

    pl_gravity = planet_gravity_for_Heliocentric_orbits(Huayno, planets_gravity.particles[0], planets_gravity.particles[1:], converter)
    gravity = Mybridge(use_threading=False)
    gravity.add_system(planetesimal_gravity, (pl_gravity,) ) ## ss works on pd
    gravity.add_system(planets_gravity, () ) # ss is self aware (evoles itself)

    orbital_parameters_for_the_planets(bodies[0:1], bodies[1:], verbose=False)
    Porb_min = orbital_period(bodies[1:].semi_major_axis.min(), bodies[0].mass)
    print "Pmin=", Porb_min
    gravity.timestep = 0.05*Porb_min

    filename = "HDdata.hdf5"
    write_set_to_file(bodies, filename, "hdf5", append_to_file=False) 

    Etot_init = planets_gravity.kinetic_energy + planets_gravity.potential_energy
    Etot = Etot_init
    dE = zero
    ddE = zero
    Mtot_init = bodies.mass.sum()
    Nenc = 0
    dEk_enc = zero
    dEp_enc = zero
    dyn_time = zero
    dt = t_end/float(n_steps)
    time = zero
    while time<t_end:
        time += dt
        bodies.age = time
        channel_from_framework_to_ssd.copy_attributes(["mass", "radius"])
        gravity.evolve_model(time)
        channel_from_ssd_to_framework.copy()
        channel_from_pdd_to_framework.copy()

#        remove_some_particles(planetesimal_gravity.orbiters_with_error, [bodies, planets_gravity.particles, planetesimal_gravity.orbiters])

        # set planetesimals baricentric
        planetesimals.position += bodies[0].position
        planetesimals.velocity += bodies[0].velocity
        Ekin = planets_gravity.kinetic_energy 
        Epot = planets_gravity.potential_energy
        ddE += (Ekin+Epot-Etot)
        Etot = Ekin + Epot
        dE = Etot_init-Etot
        Mtot = bodies.mass.sum()
        print "T=", time, t_end 
        print "M=", Mtot, "(dM[SE]=", Mtot/Mtot_init, ")",
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE/E=", dE/Etot, dE/Etot_init, "ddE/E=", ddE/Etot
        print "dE(enc)=", dEk_enc/Etot_init, dEp_enc/Etot_init
        Etot_init -= dE
        print_star(bodies[0])
        orbital_parameters_for_the_planets(bodies[0:1], bodies[1:], verbose=False)
        write_set_to_file(bodies, filename, "hdf5") 

    print_star(bodies[0])
    orbital_parameters_for_the_planets(bodies[0:1], bodies[1:], verbose=False)

def construct_single_planet_with_planetesimals(peri, apo, Ndisk):

    a = 0.5*(peri+apo)
    ecc = (apo-peri)/(apo+peri)
    print "initial orbital parameters:", a, ecc
    star_and_planet = construct_HD106906b(a, ecc)

    planetesimals = Particles(Ndisk)
    planetesimals.name = "PLN"
    planetesimals.mass = 0|units.MJupiter
    planetesimals.radius = 100|units.km

    da = (120-20)/float(Ndisk)
    a = numpy.arange(20, 120, da) 
    planetesimals.semi_major_axis = a | units.AU
    planetesimals.eccentricity = zero
    phi_rand = 2. * numpy.pi * numpy.random.rand(Ndisk)
    for i in range(len(a)):
#        planetesimals[i].position = ((a[i], 0, 0) |units.AU )
        planetesimals[i].position = (a[i] * (numpy.cos(phi_rand[i]), numpy.sin(phi_rand[i]), 0.0) |units.AU )
        vc = numpy.sqrt(constants.G*star_and_planet[0].mass/(a[i]|units.AU))
        planetesimals[i].velocity = numpy.array([-numpy.sin(phi_rand[i]), numpy.cos(phi_rand[i]), 0.0])*vc
    print planetesimals.x
    orbital_parameters_for_the_planets(star_and_planet[0:1], planetesimals, verbose=False)

    return star_and_planet, planetesimals

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", 
                      dest="n_steps", type="int", default = 1,
                      help="number of steps [%default]")
    result.add_option("-N", 
                      dest="Ndisk", type="int", default = 100,
                      help="number of disk particles [%default]")
    result.add_option("-M", unit=units.MJupiter,
                      dest="Mdisk", type="float", default = 1|units.MJupiter,
                      help="Total disk mass [%default]")
    result.add_option("-t", unit=units.yr, 
                      dest="t_end", type="float", default = 1|units.yr,
                      help="end time [%default]")
    result.add_option("-v", 
                      dest="verbose", default = False,
                      action="store_true", 
                      help="Verbosity [%default]")
    result.add_option("--peri", unit=units.AU,
                      dest="peri", type="float", default = 50|units.AU,
                      help="pericenter of planet orbit [%default]")
    result.add_option("--apo", unit=units.AU,
                      dest="apo", type="float", default = 650|units.AU,
                      help="apocenter of planet orbit [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.AU, units.yr], 
                          precision = 12, prefix = "", 
                          separator = " [", suffix = "]")

    o, arguments  = new_option_parser().parse_args()

    sun_and_planets, planetesimals = construct_single_planet_with_planetesimals(o.peri, o.apo, o.Ndisk)

    bodies = ParticlesSuperset([sun_and_planets, planetesimals])
    bodies.age = 0|units.yr

    integrate_planetary_system(sun_and_planets, planetesimals, 
                               o.t_end, o.n_steps, o.verbose)
