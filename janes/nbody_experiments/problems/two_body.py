
import math
import numpy as np
from amuse.support.data import core

from amuse.support.units import nbody_system

from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE

import sys

def particles_from_floats(m, x, v):
    stars = core.Stars(len(x))
    for i in range(len(stars)):
        star = stars[i]
        star.mass = m[i] | nbody_system.mass
        star.position = x[i] | nbody_system.length
        star.velocity = v[i] | ( nbody_system.length / nbody_system.time )
        star.radius = 0. | nbody_system.length # obligatory for Hermite() integrator
    return stars

def two_body_initial_conditions(m1, m2, a, e):
    """
    Initial conditions for a two-body problem with both bodies at the periapsis of their orbits. Assumes G=1.
    """

    m = [m1, m2]

    k1 = m2 / (m1 + m2)
    k2 = m1 / (m1 + m2)

    # http://en.wikipedia.org/wiki/Apsis
    rx_periapsis_com = (1 - e) * a
    vy_periapsis_com = math.sqrt((1 + e) * (m1 + m2) / ((1 - e) * a))
    #print vy_periapsis_com

    x = [
         [+k1 * rx_periapsis_com, 0, 0],
         [-k2 * rx_periapsis_com, 0, 0]
    ]

    v = [
         [0, +k1 * vy_periapsis_com, 0],
         [0, -k2 * vy_periapsis_com, 0]
    ]

    return particles_from_floats(m, x, v)

def two_body_orbital_period(m1, m2, a, e):
    """
    Orbital period for diagnostics with two_body_initial_conditions(m1, m2, a, e).
    """
    
    # http://www.artcompsci.org/kali/vol/two_body_problem_1/ch12.html
    m = m1 + m2
    #k1 = m2 / (m1 + m2)
    #k2 = m1 / (m1 + m2)
    #r = (1 - e) * a
    #v = math.sqrt((1 + e) * (m1 + m2) / ((1 - e) * a))
    #E = 0.5*m1*(k1*v)**2 + 0.5*m2*(k2*v)**2 - m1*m2/r
    omega3 = m / (a**2)
    omega = math.pow(omega3, 1.0/3.0)
    T = 2*math.pi / omega

    #print "Energy: %E or %E" % (E, -m1 * m2 / (2 * a))
    #print "orbital period: %s" % (T,)

    # http://en.wikipedia.org/wiki/Orbital_period (as a check)
    #print "check: %s" % (2 * math.pi * math.sqrt( math.pow(a, 3) / (m1 + m2) ), )

    return T

def binary_with_neighbours_initial_conditions():
    """
    ad-hoc binary at origin and two stars far away on nearly circular orbit
    orbital period: roughly 200
    """
    from amuse.support.data import core
    from amuse.support.units import nbody_system
    a = 10
    e = 0.9
    G = 1
    m1 = 0.5
    m2 = 0.5
    mu = G * (m1 + m2)
    h = -0.5 * mu / a
    X = 0.5*a*(1-e)
    V = 0.5*math.sqrt(-h)
    #V = 0.5 * math.sqrt( 2 * mu * ( 1 / (a * (1 - e))  - 1 / (2 * a) ) )
    x = []
    x.append([+X,  0,  0])
    x.append([-X,  0,  0])
    v = []
    v.append([ 0, -V,  0])
    v.append([ 0, +V,  0])
    x.append([+20*X,  0,  0])
    x.append([-20*X,  0,  0])
    v.append([ 0, -3*V,  0])
    v.append([ 0, +3*V,  0])
    stars = core.Stars(len(x))
    for i in range(len(stars)):
        star = stars[i]
        star.mass = 0.5 | nbody_system.mass
        star.position = x[i] | nbody_system.length
        star.velocity = v[i] | ( nbody_system.length / nbody_system.time )
        star.radius = 0. | nbody_system.length # obligatory for Hermite() integrator
    return stars

if __name__ == '__main__':
    import nbody_experiments.nbody_experiments as nbe
    import huayno2.interface as huayno_cc

    orbit_frac = 0.5
    print "orbit_frac=%s" % (orbit_frac,)
    
    # case 1: central body with an orbiting planet
    T = two_body_orbital_period(m1 = 1/20000, m2 = 10e-10, a = 4.8481e-6, e = 0.5)
    print "two_body_central period=%s" % (T, )
    r = nbe.NBodyComputation(
        initialConditions = two_body_initial_conditions(m1 = 1/20000, m2 = 10e-10, a = 4.8481e-6, e = 0.5), 
        dt = T / 100 | nbody_system.time, tFinal = orbit_frac * T | nbody_system.time, analyticSolution = None,
        ndim = 2, outfName = "two_body_central"
    )    

    # run the problem on two integrators
    r.runProblemOnIntegrator(PhiGRAPE(), label="PhiGRAPE", printProgress = False)
    r.runProblemOnIntegrator(Hermite(), label="Hermite", printProgress = False)
    hc = huayno_cc.Huayno(redirection='none')
    hc.initialize_code()
    hc.set_inttype_parameter(huayno_cc.Huayno.EVOLVE_HOLD_DKD_CC)
    hc.set_timestep_parameter(0.01)
    r.runProblemOnIntegrator(hc, label="Huayno_CC", printProgress = False)
    
    # plot a comparison in energy
    r.plotTotalEnergyRelativeError(labels = ["PhiGRAPE", "Hermite", "Huayno_CC"])

    # plot trajectories
    r.plotTrajectories(label="PhiGRAPE")
    r.plotTrajectories(label="Hermite")
    r.plotTrajectories(label="Huayno_CC")


    # case 2: equal-mass binaries
    T = two_body_orbital_period(m1 = 0.5, m2 = 0.5, a = 2.5, e = 0.5)
    print "two_body_equal period=%s" % (T, ) 
    r = nbe.NBodyComputation(
        initialConditions = two_body_initial_conditions(m1 = 0.5, m2 = 0.5, a = 2.5, e = 0.5), 
        dt = T / 100 | nbody_system.time, tFinal = orbit_frac * T | nbody_system.time, analyticSolution = None,
        ndim = 2, outfName = "two_body_equal"
    )    

    # run the problem on two integrators
    r.runProblemOnIntegrator(PhiGRAPE(), label="PhiGRAPE", printProgress = False)
    r.runProblemOnIntegrator(Hermite(), label="Hermite", printProgress = False)

    hc = huayno_cc.Huayno(redirection='none')
    hc.initialize_code()
    hc.set_inttype_parameter(huayno_cc.Huayno.EVOLVE_HOLD_DKD_CC)
    hc.set_timestep_parameter(0.0005)
    r.runProblemOnIntegrator(hc, label="Huayno_CC", printProgress = False)
    
    # plot a comparison in energy
    r.plotTotalEnergyRelativeError(labels = ["PhiGRAPE", "Hermite", "Huayno_CC"], logYScale=True)

    # plot trajectories
    r.plotTrajectories(label="PhiGRAPE")
    r.plotTrajectories(label="Hermite")
    r.plotTrajectories(label="Huayno_CC")
