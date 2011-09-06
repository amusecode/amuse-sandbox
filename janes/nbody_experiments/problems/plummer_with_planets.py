import math
import random
import numpy as np
from time import clock
from amuse.units import units
from amuse.units import nbody_system
import amuse.units.nbody_system as nbu
from amuse.ext.plummer import new_plummer_sphere, MakePlummerModel
from nbody_experiments.problems.two_body import two_body_initial_conditions, two_body_orbital_period
import nbody_experiments as nbe

from amuse import datamodel
def particles_from_floats(m, x, y, z, vx, vy, vz):
    stars = datamodel.Stars(len(x))
    for i in range(len(stars)):
        star = stars[i]
        star.mass = m[i] | nbody_system.mass
        star.position = [x[i], y[i], z[i]] | nbody_system.length
        star.velocity = [vx[i], vy[i], vz[i]] | ( nbody_system.length / nbody_system.time )
        star.radius = 0. | nbody_system.length # obligatory for Hermite() integrator
    return stars

def create_planet(m1 = 1/1000, m2 = 10e-10, a = 0.01, e = 0.0):
    #T = 2 * math.pi * math.sqrt( math.pow(a, 3) / (m1 + m2))
    #print T
    x = [(1 - e) * a, 0, 0]
    v = [0, math.sqrt((1 + e) * (m1 + m2) / ((1 - e) * a)), 0]
    return [m2, x, v]

def plummer_with_planets_initial_conditions(n=128, frac_planets=0.5, m_star = None, m_planet=10e-10, a_planet=0.001, e_planet=0.5, v=False):
    stars = new_plummer_sphere(n)
    # print average distance between stars
    if v:
        avg_star_dist = 0.0
        for stari in stars:
            for starj in stars:
                if stari != starj:
                    dvec = (starj.position - stari.position).value_in(nbu.length)
                    d = math.sqrt(dvec[0]**2 + dvec[1]**2 + dvec[2]**2)
                    avg_star_dist = avg_star_dist + d / len(stars)
        print "Average distance between stars: %f, semi-major axis: %f" % (avg_star_dist, a_planet)
        
    # normalize "integration work per planetary system", attempt #2
    if not (m_star is None):
        for star in stars:
            star.mass = m_star | nbu.mass
    planets = datamodel.Stars(int(frac_planets * n))
    for (star, planet) in zip(stars, planets):
        adj = two_body_initial_conditions(star.mass.value_in(nbu.mass), m_planet, a_planet, e_planet)
        planet.position = star.position + adj[1].position
        star.position = star.position + adj[0].position        
        planet.velocity = star.velocity + adj[1].velocity
        star.velocity = star.velocity + adj[0].velocity
        # hack, remove
        #star.mass = (1 / float(n)) | nbu.mass
        planet.mass = m_planet | nbu.mass
        planet.radius = 0. | nbu.length
        sys_orbital_period = two_body_orbital_period(star.mass.value_in(nbu.mass), planet.mass.value_in(nbu.mass), a_planet, e_planet)
        if v: print "Orbital period: %f" % (sys_orbital_period,)
    stars.add_particles(planets)
    if v: 
        print "All particles (first n/2: stars, second n/2: planets)"
        print stars
    return stars

def plummer_initial_conditions(n=128):
    particles = new_plummer_sphere(n)
    particles.scale_to_standard()
    return particles
    """
    plummer=MakePlummerModel(n)
    mass,pos,vel=plummer.new_model()
    mass=mass[0:,0]
    x=pos[0:,0]
    y=pos[0:,1]
    z=pos[0:,2]
    vx=vel[0:,0]
    vy=vel[0:,1]
    vz=vel[0:,2]
    radius=mass*0.
    tm=np.sum(mass)
    cmx=np.sum(mass*x)/tm
    cmy=np.sum(mass*y)/tm
    cmz=np.sum(mass*z)/tm
    cmvx=np.sum(mass*vx)/tm
    cmvy=np.sum(mass*vy)/tm
    cmvz=np.sum(mass*vz)/tm
    x-=cmx
    y-=cmy
    z-=cmz
    vx-=cmvx
    vy-=cmvy
    vz-=cmvz
    stars = datamodel.Stars(len(x))
    for i in range(len(stars)):
        star = stars[i]
        star.mass = mass[i] | nbody_system.mass
        star.position = [x[i], y[i], z[i]] | nbody_system.length
        star.velocity = [vx[i], vy[i], vz[i]] | ( nbody_system.length / nbody_system.time )
        star.radius = 0. | nbody_system.length # obligatory for Hermite() integrator
    return stars
    """

