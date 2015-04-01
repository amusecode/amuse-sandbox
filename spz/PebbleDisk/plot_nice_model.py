import sys, os
import numpy
import time as pytime
import platform
from numpy import random
from amuse.lab import *
from amuse.couple import bridge
from make_pebble_disk import make_pebble_disk
from amuse.units.optparse import OptionParser
from matplotlib import pyplot
from evolve_disk import *
from make_planets import *
from amuse.units import quantities

Johannes = None

def movie(stars, time, nstep):
    MEarth = 5.97219e+24 * units.kg
    pyplot.subplot(2,2,1)
    pyplot.cla()
    pyplot.subplot(2,2,2)
    pyplot.cla()
    pyplot.subplot(2,2,3)
    pyplot.cla()
    pyplot.subplot(2,2,4)
    pyplot.cla()
    m = 100*stars.mass/stars.mass.max()
    colors = ['r', 'b', 'g']
    particles = ParticlesSuperset([stars, stars[0].planets, stars[0].disk_particles])
    
    particles.move_to_center()
    print stars[0].disk_particles.mass.sum().as_quantity_in(MEarth)
    print stars[0].planets.mass / stars[0].disk_particles[0].mass 
    print stars[0].planets.mass / (1 | MEarth)
    #print stars[0]
    pyplot.subplot(2,2,1)
    pyplot.scatter(stars.x.value_in(units.AU), stars.y.value_in(units.AU), c='r', s=m)
    for star, color in zip(stars, colors):
        pyplot.scatter(star.disk_particles.x.value_in(units.AU), star.disk_particles.y.value_in(units.AU), c=color, s=1)
    
    colors = ['r', 'b', 'g','y']
    for star, color in zip(stars, colors):
        pyplot.scatter(star.planets.x.value_in(units.AU), star.planets.y.value_in(units.AU), color=colors, s=10)

    pyplot.xlabel("X [AU]")
    pyplot.ylabel("Y [AU]")
    pyplot.title("Time = "+str(time.as_quantity_in(units.yr)))

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

    pyplot.xlim(10,100)
    pyplot.ylim(0,1)
    pyplot.subplot(2,2,4)
    pyplot.ylim(0, 1)
    pyplot.semilogx()
    all_planets = Particles()
    for star in stars:
        all_planets.add_particles(star.planets)
    mscale = 0.1*all_planets.mass.min()
    for star, color in zip(stars, colors):
        pyplot.scatter(star.planets.semimajor_axis.value_in(units.AU), star.planets.eccentricity, color=colors, s=star.planets.mass/mscale, lw=0)
    pyplot.xlabel("a [AU]")
    pyplot.ylabel("e")
    pyplot.xlim(0,20)
    pyplot.ylim(0,0.2)

    pyplot.draw()
    #pyplot.savefig('xy%6.6i.png'%nstep,bbox_inches='tight')
    


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
    result.add_option("-f", 
                      dest="filename", default = "nice_model.h5",
                      help="input filename")
    return result

if __name__ == "__main__":
    o, arguments  = new_option_parser().parse_args()
    #    random.seed(seed=o.seed)
    
    if not o.filename or not os.path.exists(o.filename):
        print "Need to specifiy the name of the file to plot"
        sys.exit(1)
    stored_stars = read_set_from_file(o.filename, "amuse")
    history = list(stored_stars.iter_history())
    all = list(reversed(history))[1:]
    converter = nbody_system.nbody_to_si(stored_stars.mass.sum(), 100 | units.AU)
    Johannes = Kepler(converter)
    Johannes.initialize_code()
#    cluster_with_pebble_disks.add_system(cluster, ())
    time = 1 | units.yr
    pyplot.ion()
    nstep = 1
    for x in all[30::1]:
        stars = x.copy()
        print x.collection_attributes 
        time = x.collection_attributes.model_time
        #identify_stellar_hosts(stars)
        stars[0].disk_particles.star = stars[0]
        stars[0].planets.star = stars[0]
        #identify_stellar_hosts_for_planets(stars)
        determine_orbital_parameters_of_pebbles(stars)
        determine_orbital_parameters_of_planets(stars)
        #identify_nearest_neighbor(stars, 1000000 | units.yr)
        #pyplot.ion()
        movie(stars, time, nstep)
        nstep += 1
        #time.sleep(0.1)
    
    Johannes.stop()
    pyplot.show()
    pytime.sleep(5)
#   mencoder "mf://xy*.png" -mf fps=20 -ovc x264 -o movie.avi
