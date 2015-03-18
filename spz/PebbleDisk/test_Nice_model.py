import sys, os
import numpy
from numpy import random
from amuse.lab import *
from amuse.couple import bridge
from make_pebble_disk import make_pebble_disk
from amuse.units.optparse import OptionParser
from matplotlib import pyplot
from evolve_disk import *
from make_planets import *

Johannes = None

def movie(stars, time):
    m = 100*stars.mass/stars.mass.max()
    colors = ['r', 'b', 'g']

    pyplot.subplot(2,2,1)
    pyplot.scatter(stars.x.value_in(units.AU), stars.y.value_in(units.AU), c='r', s=m)
    for star, color in zip(stars, colors):
        pyplot.scatter(star.disk_particles.x.value_in(units.AU), star.disk_particles.y.value_in(units.AU), c=color, s=1)

    for star, color in zip(stars, colors):
        pyplot.scatter(star.planets.x.value_in(units.AU), star.planets.y.value_in(units.AU), c=color, s=10)

    pyplot.xlabel("X [AU]")
    pyplot.ylabel("Y [AU]")
    pyplot.title("Time = "+str(time))

    lim = 50
    pyplot.xlim(-lim, lim)
    pyplot.ylim(-lim, lim)

    pyplot.subplot(2,2,2)
    pyplot.scatter(stars.x.value_in(units.AU), stars.z.value_in(units.AU), c='r', s=m)
    for star, color in zip(stars, colors):
        pyplot.scatter(star.disk_particles.x.value_in(units.AU), star.disk_particles.z.value_in(units.AU), c=color, s=1)
    pyplot.xlabel("X [AU]")
    pyplot.ylabel("Z [AU]")

    pyplot.xlim(-lim, lim)
    pyplot.ylim(-lim, lim)

    pyplot.subplot(2,2,3)
    pyplot.ylim(0, 1)
    pyplot.semilogx()
    for star, color in zip(stars, colors):
        pyplot.scatter(star.disk_particles.semimajor_axis.value_in(units.AU), star.disk_particles.eccentricity, c=color, s=10, lw=0)
    pyplot.xlabel("a [AU]")
    pyplot.ylabel("e")

    pyplot.subplot(2,2,4)
    pyplot.ylim(0, 1)
    pyplot.semilogx()
    all_planets = Particles()
    for star in stars:
        all_planets.add_particles(star.planets)
    mscale = 0.1*all_planets.mass.min()
    for star, color in zip(stars, colors):
        pyplot.scatter(star.planets.semimajor_axis.value_in(units.AU), star.planets.eccentricity, c=color, s=star.planets.mass/mscale, lw=0)
    pyplot.xlabel("a [AU]")
    pyplot.ylabel("e")

    pyplot.draw()
    pyplot.subplot(2,2,1)
    pyplot.cla()
    pyplot.subplot(2,2,2)
    pyplot.cla()
    pyplot.subplot(2,2,3)
    pyplot.cla()
    pyplot.subplot(2,2,4)
    pyplot.cla()


def identify_nearest_neighbor(stars, time):

    for i, star in enumerate(stars):   
         r = (star.position - stars.position).lengths()
         r[r==(0|r.unit)] = r.max()
         istar_min = r.argmin()
         print istar_min
         star.nearest_neighbor = stars[istar_min]
         star.r_nn = r[istar_min]
         star.t_nn = time

def get_kepler_elements(star, pebble):
    pos = star.position - pebble.position
    vel = star.velocity - pebble.velocity
    Johannes.initialize_from_dyn(star.mass + pebble.mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2])
    a, e = Johannes.get_elements()
    return a, e

def determine_orbital_parameters_of_pebbles(stars):
    for star in stars:
        for pebble in star.disk_particles:
            if pebble.star is None:
                continue
            a, e = get_kepler_elements(star, pebble)
            pebble.semimajor_axis = a
            pebble.eccentricity = e

def determine_orbital_parameters_of_planets(stars):
    for star in stars:
        for planet in star.planets:
            if planet.star is None:
                continue
            a, e = get_kepler_elements(star, planet)
            planet.semimajor_axis = a
            planet.eccentricity = e

def potential_energy(pebble, star):
    r = (pebble.position-star.position).length()
    return -constants.G*pebble.mass*star.mass/r

def kinetic_energy(pebble, star):
    return 0.5*pebble.mass * (pebble.velocity-star.velocity).length()**2

def identify_host_star(pebble, stars):
    bound_to_star = None
    Eb_max = zero
    for star in stars:
        binding_energy = potential_energy(pebble, star) +  kinetic_energy(pebble, star)
        if binding_energy<zero:
            if binding_energy < Eb_max:
                Eb_max = binding_energy
                bound_to_star = star
    return bound_to_star

def identify_stellar_hosts(stars):
    for star in stars:
        star.bound_pebbles = Particles()
    for star in stars:
        for pebble in star.disk_particles:
            bound_star = identify_host_star(pebble, stars)
            pebble.star = bound_star
            if not bound_star is None:
                bound_star.bound_pebbles.add_particle(pebble)
    for star in stars:
        print "Bound pebbles:", star.mass, len(star.bound_pebbles)

def identify_stellar_hosts_for_planets(stars):
    for star in stars:
        star.bound_planets = Particles()
    for star in stars:
        for planet in star.planets:
            bound_star = identify_host_star(planet, stars)
            planet.star = bound_star
            if not bound_star is None:
                bound_star.bound_planets.add_particle(planet)
    for star in stars:
        print "Bound planets:", star.mass, len(star.bound_planets)

def remove_escapers(stars):
    for star in stars:
        pebbles = star.disk_particles
        escapers = pebbles[pebbles.position.lengths()>2*stars.position.lengths().max()]
        if len(escapers) == 0:
            return

        escapers = escapers[numpy.asarray([x is None for x in escapers.star])]
        if len(escapers)>0:
            pebbles.remove_particles(escapers)

def remove_escaping_planets(stars):
    for star in stars:
        planets = star.planets
        escapers = planets[planets.position.lengths()>2*stars.position.lengths().max()]
        if len(escapers) == 0:
            return

        escapers = escapers[numpy.asarray([x is None for x in escapers.star])]
        if len(escapers)>0:
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
                      dest="filename", default = "initial_cluster.amuse",
                      help="input filename")
    result.add_option("-R", unit=units.parsec,
                      dest="Rcluster", type="int", default = 1000|units.AU,
                      help="Cluster size [%default]")
    result.add_option("-n", 
                      dest="n", type="int", default = 1000,
                      help="Number of pebbels per star [%default]")
    return result


if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    #    random.seed(seed=o.seed)
    
    if o.filename and os.path.exists(o.filename):
        stars = read_set_from_file(o.filename, "amuse").copy()
        Rcl = o.Rcluster
        converter = nbody_system.nbody_to_si(stars.mass.sum(), Rcl)
    else:
        Ncl = o.Nstars
        masses = [1] | units.MSun
        Mcl = masses.sum()
        Rcl = o.Rcluster
        converter = nbody_system.nbody_to_si(Mcl, Rcl)
        stars = new_plummer_model(Ncl, converter)
        stars.mass = masses

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

    cluster = Hermite(converter)
    cluster.parameters.dt_param = 0.0003
#    cluster = Huayno(converter)
#    cluster.parameters.inttype_parameter = 21
#    cluster = Mercury(converter)
    from amuse.community.sakura.interface import Sakura
#    cluster = Sakura(converter)
    cluster.particles.add_particles(stars)
    channel_from_cluster = cluster.particles.new_channel_to(stars)
    tcr = Rcl/numpy.sqrt(constants.G*stars.mass.sum()/Rcl)
    print "Tcr=", tcr.in_(units.yr)
    dt = 1|units.yr
    cluster_with_pebble_disks = bridge.Bridge(timestep=dt)

    Johannes = Kepler(converter)
    Johannes.initialize_code()

    channels_from_planets = []
    for star in stars:
        cluster.particles.add_particles(star.planets)
        channels_from_planets.append(
            cluster.particles.new_channel_to(star.planets))

    for star in stars:
        pebble_gravity = advance_without_selfgravity(star.disk_particles)
        stellar_gravity = central_point_mass(star)
        cluster_with_pebble_disks.add_system(pebble_gravity, (cluster,))
        cluster_with_pebble_disks.add_system(cluster, (bridge.CalculateFieldForParticles(star.disk_particles, constants.G),))

#    cluster_with_pebble_disks.add_system(cluster, ())
    identify_stellar_hosts(stars)
    identify_stellar_hosts_for_planets(stars)
    time = 0 | units.yr
    identify_nearest_neighbor(stars, time)

    if not o.filename:
        filename = "initial_cluster.amuse"
        write_set_to_file(stars, filename, "amuse", append_to_file=False, version="2.0")
    
    pyplot.ion()
    istep = 0
    while time < 100*tcr:
        istep += 1
        time += dt
        cluster_with_pebble_disks.evolve_model(time)
        channel_from_cluster.copy()
        for channel in channels_from_planets:
            channel.copy()
        print time.in_(units.yr), cluster_with_pebble_disks.model_time.in_(units.yr)
        movie(stars, time)
        if istep%100 == 0:
            remove_escapers(stars)
            remove_escaping_planets(stars)
            identify_stellar_hosts(stars)
            identify_stellar_hosts_for_planets(stars)
            determine_orbital_parameters_of_planets(stars)
            determine_orbital_parameters_of_pebbles(stars)
            identify_nearest_neighbor(stars, time)

#            filename = "cluster{0}.amuse".format(istep)
#            write_set_to_file(stars, filename, "amuse")
    pyplot.show()
    Johannes.stop()
