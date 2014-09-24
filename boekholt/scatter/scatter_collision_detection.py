## This script runs a certain initial condition for a 
## single-binary star scattering. First with ph4 and then
## with smalln. It outputs whether there occurred a collision
## or not. 

#############################################################
## Imports
#############################################################

import sys
import numpy
import math
import time

from amuse.test.amusetest import TestWithMPI
from amuse.units import nbody_system, units, constants
from amuse.datamodel import Particles
from amuse.ic.plummer import new_plummer_model

from amuse.community.ph4.interface import ph4Interface, ph4
from amuse.community.smalln.interface import SmallN

#############################################################
## Functions
#############################################################

def get_initial_condition():
    p = Particles(3)

    p[0].mass   = 6.667e-01 | nbody_system.mass
    p[0].radius = 4.000e-03 | nbody_system.length
    p[0].x  = -1.309e+01 | nbody_system.length
    p[0].y  =  1.940e+01 | nbody_system.length
    p[0].z  = -1.163e+01 | nbody_system.length
    p[0].vx =  2.366e-01 | nbody_system.speed
    p[0].vy = -3.813e-01 | nbody_system.speed
    p[0].vz =  2.486e-01 | nbody_system.speed

    p[1].mass   = 3.333e-01 | nbody_system.mass
    p[1].radius = 1.000e-03 | nbody_system.length  
    p[1].x  = -1.506e+01 | nbody_system.length
    p[1].y  =  1.937e+01 | nbody_system.length
    p[1].z  = -1.163e+01 | nbody_system.length
    p[1].vx =  3.483e-01 | nbody_system.speed
    p[1].vy = -4.513e-01 | nbody_system.speed
    p[1].vz =  2.486e-01 | nbody_system.speed

    p[2].mass   = 5.000e-01 | nbody_system.mass
    p[2].radius = 2.000e-03 | nbody_system.length 
    p[2].x  =  2.749e+01 | nbody_system.length
    p[2].y  = -3.877e+01 | nbody_system.length
    p[2].z  =  2.325e+01 | nbody_system.length
    p[2].vx = -5.476e-01 | nbody_system.speed
    p[2].vy =  8.092e-01 | nbody_system.speed
    p[2].vz = -4.972e-01 | nbody_system.speed

    return p

def print_snapshot(t, N, t_cpu, m, pos, vel):
    print t.value_in(nbody_system.time), N, t_cpu
    i=0
    while i<N:
        print m[i].value_in(nbody_system.mass), pos[i][0].value_in(nbody_system.length), pos[i][1].value_in(nbody_system.length), pos[i][2].value_in(nbody_system.length), vel[i][0].value_in(nbody_system.speed), vel[i][1].value_in(nbody_system.speed), vel[i][2].value_in(nbody_system.speed)
        i += 1
    print " "

def get_index(particles, p):
    x = p.x.value_in(nbody_system.length)

    i=0
    while i<len(particles):
        if particles[i].x.value_in(nbody_system.length) == x:
            return i
        i += 1

    return -1

def get_ph4():
    instance = ph4()
    instance.initialize_code()
    instance.parameters.set_defaults()
    return instance

def get_smalln():
    instance = SmallN()
    instance.initialize_code()
    instance.parameters.set_defaults()
    return instance

def run(instance):
    N = 3
    t_begin = 0.0 | nbody_system.time
    t_end = 100.0 | nbody_system.time
    dt_diag = 0.125 | nbody_system.time

    particles = get_initial_condition()

    instance.particles.add_particles(particles)
    instance.commit_particles()

    kinetic_energy = instance.kinetic_energy
    potential_energy = instance.potential_energy
    E0 = kinetic_energy + potential_energy

    tcpu_tot = 0.0

    m = instance.particles.mass
    pos = instance.particles.position
    vel = instance.particles.velocity
    #print_snapshot(t_begin, N, tcpu_tot, m, pos, vel)

    sc = instance.stopping_conditions.collision_detection
    sc.enable()

    print >> sys.stderr, 'Collision detection is supported', sc.is_supported()
    print >> sys.stderr, 'Collision detection is enabled', sc.is_enabled()

    isCollision = False

    t = t_begin + dt_diag
    while t < t_end-dt_diag/2:
        tcpu0 = time.time()
 
        instance.evolve_model(t)

        tcpu = time.time()-tcpu0
        tcpu_tot += tcpu

        m = instance.particles.mass
        pos = instance.particles.position
        vel = instance.particles.velocity
        #print_snapshot(t, N, tcpu_tot, m, pos, vel)
      
        t += dt_diag

        if sc.is_set():
            print >> sys.stderr, 'collision!!'
            isCollision = True
            break

        print >> sys.stderr, t, '/', t_end, 'tcpu=', tcpu_tot

    kinetic_energy = instance.kinetic_energy
    potential_energy = instance.potential_energy
    E = kinetic_energy + potential_energy

    dE = abs((E-E0)/E0)

    print >> sys.stderr, 'tcpu=', tcpu_tot, ', dE=', dE

    instance.cleanup_code()
    instance.stop()

    return isCollision

#############################################################
## Main
#############################################################

if __name__=="__main__":

    grav_ph4 = get_ph4()
    ph4_collide = run(grav_ph4)

    grav_smalln = get_smalln()
    smalln_collide = run(grav_smalln)

    print >> sys.stderr, 'ph4 has a collision', ph4_collide
    print >> sys.stderr, 'smalln has a collision', smalln_collide





