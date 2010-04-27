import matplotlib.pyplot as plt

import sys
import unittest
import numpy
import random
import collections
import os

from amuse.support.units import nbody_system
from amuse.support.units import units

from amuse.legacy.hermite0.interface import Hermite
from amuse.legacy.bhtree.interface import BHTree

from amuse.legacy.support.core import is_mpd_running
#from amuse.test.amusetest import get_path_to_test_results

from amuse.ext.plummer import MakePlummerModel
from amuse.ext.salpeter import SalpeterIMF

import LagrangianRadii as lr

def to_com(stars) :
    com = stars.center_of_mass()
    print  "com:" , com

    vcom = stars.center_of_mass_velocity()
    print  "vcom:" , vcom
    
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
    initial_mass_function = SalpeterIMF(m_min, m_max, alpha)
    m_tot, masses = initial_mass_function.next_set(nstars)

    end_time = 10 | nbody_system.time
    dt = 0.25 | nbody_system.time 
    if not with_units :
#        dynamic_converter = Hermite.NBODY
        dynamic_converter = BHTree.NBODY
        convert_nbody = None
        masses /= m_tot.value_in(nbody_system.mass)     # scale to unit mass 
        m_tot = 1 | nbody_system.mass

    else :
        convert_nbody = nbody_system.nbody_to_si(m_tot, r_vir)
        convert_nbody.set_as_default()
        dynamic_converter = convert_nbody
        end_time = convert_nbody.to_si(end_time).as_quantity_in(units.Myr)
        dt = convert_nbody.to_si(dt).as_quantity_in(units.Myr)
        print m_tot



    stars = MakePlummerModel(nstars, convert_nbody, random_state = seed).result;
    stars.mass = masses 
    to_com(stars)
#    print "Ek=", stars.kinetic_energy()
#    print "Ep=", stars.potential_energy(gravitational_constant=1)
    
#    gravity = Hermite(dynamic_converter)
    gravity = BHTree(dynamic_converter)
    if isinstance(gravity, BHTree) :
        gravity.initialize_code()
        eps = 1/float(nstars**3) | length_unit
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
        print "Time= ", time, " E= ", Et, Ek, Ep, Ek/Ep, (Et-Et0)/Et0

    print Rl, t, r10
    plt.plot(t, r10)
    plt.plot(t, r25)
    plt.plot(t, r50)
#    axis([0, 10, 1.1*amin(r), 2*amax(r)])
    plt.axis([0, end_time.value_in(time_unit), 0, 2])
    Tlabel = 'T [' + str(time_unit) + ']'
    Rlabel = 'R_L [' + str(length_unit) + ']'
    plt.xlabel(Tlabel)
    plt.ylabel(Rlabel)
    plt.title('Lagrangian radius')
    plt.show()
