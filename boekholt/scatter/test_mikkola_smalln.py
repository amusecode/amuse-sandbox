## This script integrates an equal mass, circular binary for 
## a single period. We sample that period by 10*0.1 time-steps. 
## 
## We output the time_so_far / t_end, t_cpu, angle. 
## The angle is the angle that the orbit has progressed in
## its circular orbit. Since the time-step (for output) is
## constant, this angle should also be constant. 
##
## Issues:
## 1) Mikkola, gives two different angles, which differ by a factor 2
##    Possibly something to do with MikkolaInterface.synchronize_model
## 2) SmallN gives nan. It works with nbody_system units. 
##    Does SmallN work with physical units and converters?

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
from amuse.community.huayno.interface import Huayno
from amuse.community.hermite0.interface import Hermite
from amuse.community.sakura.interface import Sakura
from amuse.community.mikkola.interface import Mikkola
#from amuse.community.hermitepn.interface import HermitePN

#import logging
#logging.basicConfig(level=logging.DEBUG)
#logging.getLogger("code").setLevel(logging.DEBUG)

#############################################################
## Functions
#############################################################

def get_initial_condition3():
    p = Particles(3)

    p[0].mass   = 6.667e-01 | units.MSun
    p[0].radius = 4.000e-03 | units.AU
    p[0].x  = -1.309e+01 | units.AU
    p[0].y  =  1.940e+01 | units.AU
    p[0].z  = -1.163e+01 | units.AU
    p[0].vx =  2.366e-01 | units.AU / units.yr * (2.*math.pi)
    p[0].vy = -3.813e-01 | units.AU / units.yr * (2.*math.pi)
    p[0].vz =  2.486e-01 | units.AU / units.yr * (2.*math.pi)

    p[1].mass   = 3.333e-01 | units.MSun
    p[1].radius = 1.000e-03 | units.AU
    p[1].x  = -1.506e+01 | units.AU
    p[1].y  =  1.937e+01 | units.AU
    p[1].z  = -1.163e+01 | units.AU
    p[1].vx =  3.483e-01 | units.AU / units.yr * (2.*math.pi)
    p[1].vy = -4.513e-01 | units.AU / units.yr * (2.*math.pi)
    p[1].vz =  2.486e-01 | units.AU / units.yr * (2.*math.pi)

    p[2].mass   = 5.000e-01 | units.MSun
    p[2].radius = 2.000e-03 | units.AU
    p[2].x  =  2.749e+01 | units.AU
    p[2].y  = -3.877e+01 | units.AU
    p[2].z  =  2.325e+01 | units.AU
    p[2].vx = -5.476e-01 | units.AU / units.yr * (2.*math.pi)
    p[2].vy =  8.092e-01 | units.AU / units.yr * (2.*math.pi)
    p[2].vz = -4.972e-01 | units.AU / units.yr * (2.*math.pi)

    return p

def get_initial_condition2():
    p = Particles(2)

    p[0].mass   = 0.5 | units.MSun
    p[0].radius = 1e-3 | units.AU
    p[0].x  = 0.5 | units.AU
    p[0].y  = 0. | units.AU
    p[0].z  = 0. | units.AU
    p[0].vx = 0. | units.AU / units.yr * (2.*math.pi)
    p[0].vy = 0.5 | units.AU / units.yr * (2.*math.pi)
    p[0].vz = 0. | units.AU / units.yr * (2.*math.pi)

    p[1].mass   = 0.5 | units.MSun
    p[1].radius = 1e-3 | units.AU
    p[1].x  = -0.5 | units.AU
    p[1].y  =  0. | units.AU
    p[1].z  =  0. | units.AU
    p[1].vx =  0. | units.AU / units.yr * (2.*math.pi)
    p[1].vy = -0.5 | units.AU / units.yr * (2.*math.pi)
    p[1].vz =  0. | units.AU / units.yr * (2.*math.pi)

    return p

def print_snapshot(t, N, t_cpu, m, pos, vel):
    print t.number, N, t_cpu
    i=0
    while i<N:
        print m[i].value_in(units.MSun), pos[i][0].value_in(units.AU), pos[i][1].value_in(units.AU), pos[i][2].value_in(units.AU), vel[i][0].value_in(units.AU/units.yr)/(2.*math.pi), vel[i][1].value_in(units.AU/units.yr)/(2.*math.pi), vel[i][2].value_in(units.AU/units.yr)/(2.*math.pi)
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

def get_ph4(converter):
    instance = ph4(converter)
    instance.initialize_code()
    instance.parameters.set_defaults()
    return instance

def get_smalln(converter):
    instance = SmallN(converter)
    instance.initialize_code()
    instance.parameters.set_defaults()
    return instance

def get_huayno(converter):
    instance = Huayno(converter)
    instance.initialize_code()
    instance.parameters.set_defaults()
    return instance

def get_hermite(converter):
    instance = Hermite(converter)
    instance.initialize_code()
    instance.parameters.set_defaults()
    return instance

def get_sakura(converter):
    instance = Sakura(converter)
    instance.initialize_code()
    instance.parameters.set_defaults()
    instance.parameters.timestep = 0.1 | units.yr
    return instance

def get_mikkola(converter):
    instance = Mikkola(converter)
    instance.initialize_code()
    instance.parameters.set_defaults()
    instance.parameters.lightspeed = constants.c
    #instance.parameters.timestep = 0.01 | nbody_system.time
    #print >> sys.stderr, 'c/v=', instance.parameters.lightspeed

    print >> sys.stderr, instance.parameters

    return instance

#def get_hermitepn(converter):
#    instance = HermitePN(converter)
#    instance.initialize_code()
#    instance.parameters.set_defaults()
#    return instance

def check_for_collision_pericenter(mass, radius, pos, vel):
    N = len(mass)
    i=0
    while i<N-1:
        j=i+1
        while j<N:
            mu  = (mass[i]+mass[j]).value_in(units.MSun)
            dx  = (pos[i][0]-pos[j][0]).value_in(units.AU)
            dy  = (pos[i][1]-pos[j][1]).value_in(units.AU)
            dz  = (pos[i][2]-pos[j][2]).value_in(units.AU)
            dvx = (vel[i][0]-vel[j][0]).value_in(units.AU/units.yr)/(2.*math.pi)
            dvy = (vel[i][1]-vel[j][1]).value_in(units.AU/units.yr)/(2.*math.pi)
            dvz = (vel[i][2]-vel[j][2]).value_in(units.AU/units.yr)/(2.*math.pi)

            dr2 = dx**2 + dy**2 + dz**2
            dr = math.sqrt(dr2)

            dv2 = dvx**2 + dvy**2 + dvz**2

            rdotv = dx*dvx + dy*dvy + dz*dvz

            e = 0.5*dv2 - mu/dr

            lx = dy*dvz - dz*dvy
            ly = dz*dvx - dx*dvz
            lz = dx*dvy - dy*dvx
            l2 = lx**2 + ly**2 + lz**2

            if e < 0.:
                D = mu**2 + 2.*e*l2
                rp = (-mu+math.sqrt(D)) / (2.0*e)
            elif e == 0.:
                rp = l2/(2.*mu)
            else:
                D = mu**2 + 2.*e*l2
                rp = (-mu+math.sqrt(D)) / (2.0*e)

            if rp < (radius[i]+radius[j]).value_in(units.AU):
                return True, rdotv

            j += 1
        i += 1
    return False, rdotv

def get_angle(r, r0):
    x = r[1][0]-r[0][0]
    y = r[1][1]-r[0][1]
    z = r[1][2]-r[0][2]
    
    x0 = r0[1][0]-r0[0][0]
    y0 = r0[1][1]-r0[0][1]
    z0 = r0[1][2]-r0[0][2]

    dot = (x*x0 + y*y0 + z*z0).number

    r2 = (x**2 + y**2 + z**2).number
    r02 = (x0**2 + y0**2 + z0**2).number
 
    cos_angle = dot / ( math.sqrt(r2)*math.sqrt(r02) )

    return math.acos(cos_angle)

def run(instance):
    t_begin =  0.0 | units.yr
    t_end   =  1.0 | units.yr
    dt_diag =  0.1 | units.yr

    N = 2
    particles = get_initial_condition2()

    instance.particles.add_particles(particles)
    instance.commit_particles()

    kinetic_energy = instance.kinetic_energy
    potential_energy = instance.potential_energy
    E0 = kinetic_energy + potential_energy

    tcpu_tot = 0.0

    m = instance.particles.mass
    r = instance.particles.radius
    pos = instance.particles.position
    vel = instance.particles.velocity
    #print_snapshot(t_begin, N, tcpu_tot, m, pos, vel)

    try: 
        sc = instance.stopping_conditions.collision_detection
        sc.enable()
        sc_is_enabled = True
        print >> sys.stderr, 'Collision detection is supported', sc.is_supported()
        print >> sys.stderr, 'Collision detection is enabled', sc.is_enabled()
    except:
        sc_is_enabled = False
        print >> sys.stderr, 'Collision detection is not supported'
        print >> sys.stderr, 'Collision detection is not enabled'

    print >> sys.stderr, sc_is_enabled
    isCollision = False
    rdotv_prev = 0.

    t = t_begin + dt_diag
    while t < t_end-dt_diag/2:
        tcpu0 = time.time()

        r_prev = instance.particles.position

        try: 
            instance.evolve_model(t)
        except:
            m = instance.particles.mass
            r = instance.particles.radius
            pos = instance.particles.position
            vel = instance.particles.velocity
            #print_snapshot(t, N, tcpu_tot, m, pos, vel)
            print >> sys.stderr, 'collision by except evolve!!'
            isCollision = True
            break

        tcpu = time.time()-tcpu0
        tcpu_tot += tcpu

        r_now = instance.particles.position
        angle = get_angle(r_now, r_prev)

        m = instance.particles.mass
        r = instance.particles.radius
        pos = instance.particles.position
        vel = instance.particles.velocity
        #print_snapshot(t, N, tcpu_tot, m, pos, vel)
        #print >> sys.stderr, t, m, r
        t += dt_diag

        if sc_is_enabled:
            #print >> sys.stderr, 'checking sc is set...'
            if sc.is_set():
                print >> sys.stderr, 'collision by sc!!'
                isCollision = True
                break
        else:
            #print >> sys.stderr, 'checking manually...'
            coll, rdotv = check_for_collision_pericenter(m, r, pos, vel)
            #print >> sys.stderr, coll, rdotv, rdotv_prev
            if coll:
                if rdotv_prev < 0 and rdotv > 0:
                    print >> sys.stderr, 'collision by check!!'
                    isCollision = True
                    break
                else:
                    rdotv_prev = rdotv
                    #print >> sys.stderr, rdotv_prev, rdotv
            else:  
                rdotv_prev = rdotv

        #print >> sys.stderr, instance.t
        print >> sys.stderr, t, '/', t_end, 'tcpu=', tcpu_tot, 'angle=', angle

    m = instance.particles.mass
    r = instance.particles.radius
    pos = instance.particles.position
    vel = instance.particles.velocity
    #print_snapshot(t, N, tcpu_tot, m, pos, vel)

    kinetic_energy = instance.kinetic_energy
    potential_energy = instance.potential_energy
    E = kinetic_energy + potential_energy

    dE = abs((E-E0)/E0)

    print >> sys.stderr, 'tcpu=', tcpu_tot, ', dE=', dE

    instance.cleanup_code()
    instance.stop()

    return t, isCollision, tcpu_tot, dE

#############################################################
## Main
#############################################################

if __name__=="__main__":

    converter = nbody_system.nbody_to_si(1|units.MSun,1|units.AU)

    print >> sys.stderr, ' '
    print >> sys.stderr, ' '

    print >> sys.stderr, 'ph4 ...'
    print >> sys.stderr, ' '
    grav_ph4 = get_ph4(converter)
    t_ph4, ph4_collide, tcpu_ph4, dE_ph4 = run(grav_ph4)

    print >> sys.stderr, 'smalln ...'
    print >> sys.stderr, ' '
    grav_smalln = get_smalln(converter)
    t_smalln, smalln_collide, tcpu_smalln, dE_smalln = run(grav_smalln)

    #print >> sys.stderr, 'huayno ...'
    #grav_huayno = get_huayno(converter)
    #t_huayno, huayno_collide, tcpu_huayno, dE_huayno = run(grav_huayno)

    #print >> sys.stderr, 'hermite ...'
    #grav_hermite = get_hermite(converter)
    #t_hermite, hermite_collide, tcpu_hermite, dE_hermite = run(grav_hermite)

    print >> sys.stderr, 'sakura ...'
    print >> sys.stderr, ' '
    grav_sakura = get_sakura(converter)
    t_sakura, sakura_collide, tcpu_sakura, dE_sakura = run(grav_sakura)

    print >> sys.stderr, 'mikkola ...'
    print >> sys.stderr, ' '
    grav_mikkola = get_mikkola(converter)
    t_mikkola, mikkola_collide, tcpu_mikkola, dE_mikkola = run(grav_mikkola)

    #print >> sys.stderr, 'hermitepn ...'
    #grav_hermitepn = get_hermitepn(converter)
    #t_hermitepn, hermitepn_collide, tcpu_hermitepn, dE_hermitepn = run(grav_hermitepn)

    print >> sys.stderr, 'ph4 has a collision', ph4_collide, 'at', t_ph4, 'which took', tcpu_ph4, 'dE=', dE_ph4
    print >> sys.stderr, 'smalln has a collision', smalln_collide, 'at', t_smalln, 'which took', tcpu_smalln, 'dE=', dE_smalln
    #print >> sys.stderr, 'huayno has a collision', huayno_collide, 'at', t_huayno, 'which took', tcpu_huayno, 'dE=', dE_huayno
    #print >> sys.stderr, 'hermite has a collision', hermite_collide, 'at', t_hermite, 'which took', tcpu_hermite, 'dE=', dE_hermite
    print >> sys.stderr, 'sakura has a collision', sakura_collide, 'at', t_sakura, 'which took', tcpu_sakura, 'dE=', dE_sakura
    print >> sys.stderr, 'mikkola has a collision', mikkola_collide, 'at', t_mikkola, 'which took', tcpu_mikkola, 'dE=', dE_mikkola
    #print >> sys.stderr, 'hermitepn has a collision', hermitepn_collide, 'at', t_hermitepn, 'which took', tcpu_hermitepn, 'dE=', dE_hermitepn

    print >> sys.stderr, ' '
    print >> sys.stderr, ' '








