run_with_Hermite = False

import matplotlib.pyplot as plt

import sys
import unittest
import numpy
import random
import collections
import os


from amuse.units import nbody_system
from amuse.units import units
from amuse.community.hermite0.interface import Hermite
from amuse.community.bhtree.interface import BHTree

import LagrangianRadii as lr
from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution
"""
   HermiteOrBHTree.py
   Example code for running a simple N-body system either using either
   the hermite0 or the BHTree code for integrating the equations of
   motion.
   Written in Spring 2010 by Simon Portegies Zwart <spz@strw.leidenuniv.nl>
"""

def scale(stars) :
    # function scales positions, such that Potential Energy = -1/2
    # and scales velocities, such that Kinetic Energy = 1/4
    Ek = stars.kinetic_energy()
    Ep = stars.potential_energy(G=nbody_system.G)
    Et = Ek + Ep
    print "Time= 0, E= ", Et, Ek, Ep, " Q= ", Ek/Ep


    stars.position = stars.position * (-2.0*Ep.number)
    stars.velocity = stars.velocity / numpy.sqrt(4.0*Ek.number)

def to_com(stars) :
    com = stars.center_of_mass()
    print  "com:", com

    vcom = stars.center_of_mass_velocity()
    print  "vcom:", vcom
    
    stars.position = stars.position - com
    stars.velocity = stars.velocity - vcom
    
if __name__ == '__main__':

    assert is_mpd_running()
    seed = None
#   seed = numpy.random.RandomState([1,1,1])

    nstars = 128
    if len(sys.argv)>1 :
        nstars = int(sys.argv[1])
    with_units = len(sys.argv) > 2

    if not with_units :
        mass_unit = nbody_system.mass
        length_unit = nbody_system.length
        time_unit = nbody_system.time
    else :
        mass_unit = units.MSun
        length_unit = units.parsec
        time_unit = units.Myr

    m_min = 0.1 | mass_unit
    m_max = 100 | mass_unit
    alpha = -2.35

    r_vir = 1 | length_unit
    masses = new_salpeter_mass_distribution(nstars, m_min, m_max, alpha)
    m_tot = masses.sum()

    end_time = 1 | nbody_system.time
    dt = 0.25 | nbody_system.time 
    if not with_units :
        convert_nbody = None
        masses /= m_tot.value_in(nbody_system.mass)     # scale to unit mass 
        m_tot = 1 | nbody_system.mass

    else :
        convert_nbody = nbody_system.nbody_to_si(m_tot, r_vir)
        convert_nbody.set_as_default()
        end_time = convert_nbody.to_si(end_time).as_quantity_in(units.Myr)
        dt = convert_nbody.to_si(dt).as_quantity_in(units.Myr)
        print m_tot

    stars = new_plummer_model(nstars, convert_nbody, random_state = seed);
    stars.mass = masses 
    to_com(stars)
    scale(stars)

    if run_with_Hermite :
        gravity = Hermite()
    else :
        gravity = BHTree()
    gravity.initialize_code()
    # Set the time step constant
    gravity.legacy_interface.dt_param = 0.03
    if isinstance(gravity, BHTree) :
        gravity.parameters.timestep = 0.004 | time_unit
        eps = 1.0 / 20.0 / (nstars**0.33333) | length_unit 
    else :
        eps = 0 | length_unit
    gravity.parameters.epsilon_squared = eps**2 

    stars.radius = 0.0 | length_unit

    gravity.particles.add_particles(stars)
    gravity.commit_particles()
    from_model_to_gravity = stars.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(stars)
    
    time = 0.0 | time_unit
    Ek0 = gravity.kinetic_energy
    Ep0 = gravity.potential_energy
    Et0 = Ek0 + Ep0
    print "Time= ", time, " E= ", Et0, Ek0, Ep0, " Q= ", Ek0/Ep0

    print "evolving the model until t = " + str(end_time)
    t = []
    r10 = []
    r25 = []
    r50 = []
    r75 = []
    while time < end_time:
        time += dt
        
        gravity.evolve_model(time)
        from_gravity_to_model.copy()

        Ek = gravity.kinetic_energy
        Ep = gravity.potential_energy
        Et = Ek + Ep

        Rl = lr.LagrangianRadii(stars)
        t.append(time.value_in(time_unit))
        r10.append(Rl[4].value_in(length_unit))
        r25.append(Rl[5].value_in(length_unit))
        r50.append(Rl[6].value_in(length_unit))
        r75.append(Rl[7].value_in(length_unit))
        print "Time = ", time, \
            "\n   E = ", Et, "\n  Ek = ", Ek, "\n  Ep = ", Ep, \
            "\n   Q = ", Ek/Ep, (Et-Et0)/Et0

#    print Rl, t, r10
    plt.plot(t, r10)
    plt.plot(t, r25)
    plt.plot(t, r50)
    plt.plot(t, r75)
#    axis([0, 10, 1.1*amin(r), 2*amax(r)])
    plt.axis([0, end_time.value_in(time_unit), 0, 2])
    Tlabel = 'T [' + str(time_unit) + ']'
    Rlabel = 'R_L [' + str(length_unit) + ']'
    plt.xlabel(Tlabel)
    plt.ylabel(Rlabel)
    plt.title('Lagrangian radius')
    plt.show()
