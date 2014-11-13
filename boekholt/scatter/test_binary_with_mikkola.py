## This script integrates an equal mass, circular binary for 
## a single period. We sample that period by 10*0.1 time-steps. 
## 
## We output whether the true anomaly is integrated correctly
## with equal angles in equal time-intervals.  
## The angle is the angle that the orbit has progressed in
## its circular orbit. Since the time-step is
## constant, this angle should also be constant. 
##
## Issues:
## 1) Mikkola, gives two different angles, which differ by a factor 2
##    Possibly something to do with MikkolaInterface.synchronize_model

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

def get_initial_condition():
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

def get_ph4(converter):
    instance = ph4(converter)
    instance.initialize_code()
    instance.parameters.set_defaults()
    return instance

def get_mikkola(converter):
    instance = Mikkola(converter)
    instance.initialize_code()
    instance.parameters.set_defaults()
    instance.parameters.lightspeed = constants.c
    #instance.parameters.timestep = 0.01 | nbody_system.time
    #print >> sys.stderr, 'c/v=', instance.parameters.lightspeed

    return instance

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
    particles = get_initial_condition()

    instance.particles.add_particles(particles)
    instance.commit_particles()

    angles = []

    t = t_begin + dt_diag
    while t < t_end-dt_diag/2:
        r_prev = instance.particles.position

        instance.evolve_model(t)

        r_now = instance.particles.position
        angle = get_angle(r_now, r_prev)
        angles.append(angle)

        t += dt_diag

    print >> sys.stderr, 'correct evolution =',
    #print >> sys.stderr, 'tcpu=', tcpu_tot, ', dE=', dE, ', correct evolution =',

    instance.cleanup_code()
    instance.stop()

    angle_mean = 1.*sum(angles)/len(angles)
    angle_var  = [(a-angle_mean)**2 for a in angles]
    angle_var  = sum(angle_var)/(len(angles)-1.)
    angle_sd   = math.sqrt(angle_var)

    # Allow for a certain spread in what should be constant angles
    f = angle_sd / angle_mean
    if f < 0.1:
        return True
    else:
        return False

#############################################################
## Main
#############################################################

if __name__=="__main__":

    converter = nbody_system.nbody_to_si(1|units.MSun,1|units.AU)

    print >> sys.stderr, ' '

    grav_ph4 = get_ph4(converter)
    print >> sys.stderr, 'ph4 ...', run(grav_ph4)

    grav_mikkola = get_mikkola(converter)
    print >> sys.stderr, 'mikkola ...', run(grav_mikkola)

    print >> sys.stderr, ' '








