"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
import numpy
from  time import clock
from amuse.lab import *
from Kepler.interface import Kepler
from Kepler_orbiters import Orbiters
from Kepler_orbiters import orbital_period, orbital_parameters_for_the_planets, construct_orbital_orientation, reconstruct_posvel

from amuse.ext.stellar_tidal_evolution import TidalEvolution
from amuse.ext.stellar_gyration_radius import calculate_gyration_radius
from amuse.couple import bridge
from amuse.units import quantities

from starplanet import get_kepler_elements, orbital_period, construct_HD106906b

He_WD = 10 | units.stellar_type

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

"""
def construct_planetsimal_disk(Rinner, Router, Ndisk, Mdisk):
    planetesimals = Particles(Ndisk)
    planetesimals.semi_major_axis = Rinner + (Router-Rinner)*random.random()
    planetesimals.eccentricity = 0
    planetesimals.mass = Mdisk/float(Ndisk)
    return planetesimals
"""

def construct_planetesimal_disk(Mstar, Rinner, Router, Ndisk, Mdisk):

    planetesimals = Particles(Ndisk)
    planetesimals.mass = Mdisk/float(Ndisk)
    planetesimals.name = "PLN"
    planetesimals.radius = 100|units.km
    for pi in planetesimals:
        a = Rinner +  (Router-Rinner) * numpy.random.random() 
        pi.position = (1, 0, 0)*a
        vc = numpy.sqrt(constants.G*Mstar/a)
        pi.velocity = (0, 1, 0)*vc
    return planetesimals

def construct_planetarysystem(a, ecc, Mdisk, Ndisk, Rinner, Router):

    sun_and_planet = construct_HD106906b(a, ecc)
    planetesimals = construct_planetesimal_disk(sun_and_planet[0].mass, Rinner, Router, Ndisk, Mdisk)

    return sun_and_planet, planetesimals

def unbound_particles(bodies):
    lost_particles = Particles()
    for bi in bodies:
        a = bi.semi_major_axis
        e = bi.eccentricity
        if e<0 or e>=1:
            print "Escaper:", bi.name, "a=", a, "e=", e
            lost_particles.add_particle(bi)
    return lost_particles

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
    
def print_solar_system(bodies):
    from evolve_star import print_star
    print_star(bodies[0])
    for bi in bodies[1:8]:
        print "Planet: ", bi.age, bi.name, "a=", bi.semi_major_axis, "e=", bi.eccentricity
    for bi in bodies[9:]:
        print "Planetesimal: ", bi.age, bi.name, "a=", bi.semi_major_axis, "e=", bi.eccentricity

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
    else:
        print "No particles removed"

def starplanetwithdisk(peri, apo, t_end, n_steps, Rinner, Router, Mdisk, Ndisk):

    restartfile = "RestartHD.hdf5"
    a = 0.5*(peri+apo)
    ecc = (apo-peri)/(apo+peri)
    print "generate new initial conditions"
    sun_and_planets, planetesimals = construct_planetarysystem(a, ecc, Mdisk, Ndisk, Rinner, Router)
    bodies = ParticlesSuperset([sun_and_planets, planetesimals])
    t_start = 13|units.Myr
    bodies.age = t_start
    bodies.semi_major_axis = 0|units.RSun
    bodies.eccentricity = 0.0
    bodies.add_vector_attribute('longitudinal_unit_vector', ('lx','ly','lz')) 
    sun_and_planets.add_vector_attribute('longitudinal_unit_vector', ('lx','ly','lz')) 
    bodies.add_vector_attribute('transverse_unit_vector', ('tx','ty','tz')) 
    sun_and_planets.add_vector_attribute('transverse_unit_vector', ('tx','ty','tz')) 

    sun = sun_and_planets[0:1]
    planets = sun_and_planets - sun
    planetesimals = bodies - sun_and_planets
    orbital_parameters_for_the_planets(sun, bodies-sun)
    construct_orbital_orientation(sun, planets)

    converter=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)

    # add specific stellar parameters
    bodies[0].gyration_radius_sq = 0.2
    bodies[0].Omega = 2.6e-6|units.s**-1

#    stellar = MESA()
#    stellar.parameters.RGB_wind_scheme=1 
#    stellar.parameters.AGB_wind_scheme=1 
    stellar = SeBa()
#    stellar.parameters.RGB_wind_scheme=1 
#    stellar.parameters.AGB_wind_scheme=1 
    stellar.parameters.metallicity = 0.02
    stellar.particles.add_particles(bodies[0:1])

    channel_from_se_to_framework = stellar.particles.new_channel_to(bodies)

    print "evolve star to: ", t_start
    if t_start>zero:
        stellar.evolve_model(t_start)
    #sun.gyration_radius_sq = calculate_gyration_radius(stellar.particles)
    sun.Omega = sun.Omega * (sun.radius/stellar.particles[0].radius)**2
    channel_from_se_to_framework.copy_attributes(["mass", "age", "radius", "stellar_type", "luminosity", "temperature", "stellar_type"])

    gravity = Huayno(converter)
    gravity.particles.add_particles(sun_and_planets)
    channel_from_gravity_to_framework = gravity.particles.new_channel_to(bodies)
    channel_from_framework_to_gravity = bodies.new_channel_to(gravity.particles)

    orbital_parameters_for_the_planets(sun, planets)

    orbiter = Orbiters(converter)
    orbiter.setup_model(planetesimals, sun)
    orbiter.particles = orbiter.orbiters
    orbiter.dmdt = 0 | units.MSun/units.yr
    orbiter.dt = t_start

    channel_from_orb_to_framework = orbiter.orbiters.new_channel_to(bodies)
    channel_from_framework_to_orb = bodies.new_channel_to(orbiter.orbiters)
    channel_from_cp_to_framework = orbiter.central_particle.new_channel_to(bodies)
    channel_from_framework_to_cp = bodies.new_channel_to(orbiter.central_particle)

    channel_from_orb_to_framework.copy_attributes(["semi_major_axis", "eccentricity"])

    tidal = TidalEvolution(sun)
    tidal.add_particles(bodies-sun)
    channel_from_tidal_p = tidal.orbiters.new_channel_to(bodies)
    channel_to_tidal_p = bodies.new_channel_to(tidal.orbiters)
    channel_from_tidal_s = tidal.central_particle.new_channel_to(bodies)
    channel_to_tidal_s = bodies.new_channel_to(tidal.central_particle)

    lost_particles = unbound_particles(bodies)
    print "N unbound:", len(lost_particles)
    remove_some_particles(lost_particles, [bodies, gravity.particles, orbiter.orbiters, tidal.orbiters])

    channel_to_tidal_s.copy_attributes(["mass", "radius", "gyration_radius_sq", "Omega"])

    st_prev = sun.stellar_type
    istp = st_prev.value_in(units.stellar_type)
    if not restartfile:
        filename = 'SolarSystem_%(i)04d.hdf5' %{'i':istp}
        write_set_to_file(bodies, filename, "hdf5", append_to_file=False) 

    # First time write output file
    write_set_to_file(bodies, 'SolarSystem_data.hdf5', "hdf5", append_to_file=False) 

    time = sun.age
    print "evolved to: T=", time, t_start, "M=", sun.mass
    Mtot_init = bodies[0].mass
    while sun.stellar_type<He_WD:
        write_set_to_file(bodies, "RestartSolarSystem.hdf5", "hdf5", 
                          append_to_file=False) 

        t_prev = sun.age
        m_prev = sun.mass
        r_prev = sun.radius
        #stellar.evolve_model()

        #sun.gyration_radius_sq = calculate_gyration_radius(stellar.particles)
        channel_from_se_to_framework.copy_attributes(["age", "mass", "radius", "luminosity", "temperature", "stellar_type"])
        sun.Omega = sun.Omega * (r_prev/sun.radius)**2
        dm = m_prev-sun.mass
        channel_to_tidal_s.copy_attributes(["mass", "radius", "gyration_radius_sq", "Omega", "luminosity", "temperature", "stellar_type"])

        orbital_parameters_for_the_planets(sun, planets)
        construct_orbital_orientation(sun, planets)

        channel_to_tidal_p.copy_attributes(["semi_major_axis", "eccentricity"])
        channel_to_tidal_s.copy_attributes(["Omega"])

        tidal.evolve_model(time)

        channel_from_tidal_p.copy_attributes(["semi_major_axis", "eccentricity"])
        channel_from_tidal_s.copy_attributes(["Omega"])

        reconstruct_posvel(sun, planets)
        channel_from_gravity_to_framework.copy_attributes(["x", "y", "z",
                                                           "vx", "vy", "vz"])

        remove_some_particles(tidal.all_merged_orbiters, [bodies, gravity.particles, orbiter.orbiters, tidal.orbiters])
        remove_some_particles(tidal.orbiters_with_error, [bodies, gravity.particles, orbiter.orbiters, tidal.orbiters])

        if len(planetesimals)==0:
            print "No planetesimals left: exit"
            exit()

        if len(planets)==0:
            print "No planetesimals left: exit"
            exit()

        time = bodies[0].age
        dm = bodies[0].mass-m_prev
        dmdt = 0 | units.MSun/units.yr

        m_post = bodies[0].mass
        bodies[0].mass = m_prev
        channel_from_framework_to_cp.copy_attributes(["mass"])
        dt = time-t_prev

        print "dmdt=", dmdt, time, dt

        verbose = True
        trun_start = clock() | units.s
        one_step_for_planetary_system(sun, planets, planetesimals, 
                                      gravity, orbiter,  converter,
                                      time, time+dt, n_steps, verbose)
        trun_stop = clock() | units.s
        print "Performace data: N=", len(sun), len(planets), "n=", len(planetesimals), "t_end=", t_end, time, "dt for step=", trun_stop - trun_start
        time += dt

        remove_some_particles(gravity.orbiters_with_error, [bodies, gravity.particles, orbiter.orbiters, tidal.orbiters])
               
        channel_from_cp_to_framework.copy_attributes(["mass", "x", "y", "z","vx", "vy", "vz"])
        channel_from_orb_to_framework.copy_attributes(["x", "y", "z","vx", "vy", "vz", "semi_major_axis", "eccentricity"])

        bodies.position -= sun.position
        bodies.velocity -= sun.velocity
        channel_from_framework_to_gravity.copy_attributes(["x", "y", "z",
                                                           "vx", "vy", "vz"])
        orbital_parameters_for_the_planets(sun, planets)
        #print_solar_system(bodies)
 
        bodies[1:].age = gravity.time

        # dump system
        stp = bodies[0].stellar_type
        if not st_prev == stp:
            print "STP:", stp, st_prev
            st_prev = stp
            istp = stp.value_in(units.stellar_type)
            filename = 'SolarSystem_%(i)04d.hdf5' %{'i':istp}
            write_set_to_file(bodies, filename, "hdf5", append_to_file=False)

        Mtot = bodies.mass.sum()
        print "T=", time
        print "N=", len(bodies)
        print "M=", Mtot, "(dM[SE]=", Mtot/Mtot_init, ") MSun=", bodies[0].mass
        print "Masses:", bodies[0].mass, m_prev, m_post 
        #print_solar_system(bodies)
        write_set_to_file(bodies, 'SolarSystem_data.hdf5', "hdf5", append_to_file=True) 

    gravity.stop()
    stellar.stop()

def one_step_for_planetary_system(sun, planets, planetesimals, 
                                  planets_gravity, planetesimal_gravity, 
                                  converter, 
                                  t_start, t_end, n_steps, verbose):

    sun_and_planets = ParticlesSuperset([sun, planets])
    bodies = ParticlesSuperset([sun, planets, planetesimals])

    channel_from_ssd_to_framework = planets_gravity.particles.new_channel_to(sun_and_planets)
    channel_from_framework_to_ssd = sun_and_planets.new_channel_to(planets_gravity.particles)
    channel_from_pdd_to_framework = planetesimal_gravity.particles.new_channel_to(planetesimals)
    channel_from_framework_to_pdd = planetesimals.new_channel_to(planetesimal_gravity.particles)

    pl_gravity = planet_gravity_for_Heliocentric_orbits(Huayno, planets_gravity.particles[0], planets_gravity.particles[1:], converter)

    gravity = Mybridge(use_threading=False)
    Porb_min = orbital_period(planets.semi_major_axis.min(), sun.mass)
    print "Pmin=", Porb_min
    gravity.timestep = 0.05*Porb_min
    gravity.add_system(planetesimal_gravity, (pl_gravity,) ) ## ss works on pd
    gravity.add_system(planets_gravity, () ) # ss is self aware (evoles itself)

    channel_from_ssd_to_framework.copy_attributes(["x", "y", "z",
                                                   "vx", "vy", "vz"])
    channel_from_pdd_to_framework.copy_attributes(["x", "y", "z",
                                                   "vx", "vy", "vz", 
                                                   "semi_major_axis", 
                                                   "eccentricity"])

    Etot_init = planets_gravity.kinetic_energy + planets_gravity.potential_energy
    Etot = Etot_init
    dE = zero
    ddE = zero
    Mtot_init = bodies.mass.sum()
    Nenc = 0
    dEk_enc = zero
    dEp_enc = zero
    dt = (t_end-t_start)/float(n_steps)
    time = t_start
    while time<t_end:
        time += dt
        gravity.evolve_model(time)
        bodies.age = gravity.time

        channel_from_ssd_to_framework.copy_attributes(["x", "y", "z",
                                                       "vx", "vy", "vz"])
        channel_from_pdd_to_framework.copy_attributes(["x", "y", "z",
                                                       "vx", "vy", "vz", 
                                                       "semi_major_axis", 
                                                       "eccentricity"])
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

        orbital_parameters_for_the_planets(sun, planets, verbose=verbose)
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="int", default = 10,
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
                      dest="t_end", type="float", default = 100|units.yr,
                      help="end time of the simulation [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.AU, units.yr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")
    o, arguments  = new_option_parser().parse_args()
    starplanetwithdisk(**o.__dict__)

