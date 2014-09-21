"""
testing orbital elements in Kepler
"""

from amuse.units import units, constants, nbody_system
from amuse.community.kepler.interface import Kepler
from amuse.datamodel import Particles

import numpy

sun_and_stone = Particles(2)
sun_and_stone[0].position = [-5.40085336308e+13, -5.92288255636e+13, 0.] | units.m
sun_and_stone[0].velocity = [-68.4281023588, -90.213640012, 0.] | (units.m/units.s)
sun_and_stone[0].mass = 1.98892e+30 | units.kg
sun_and_stone[1].position = [2.31421952844e+17, 1.72742642121e+17, 0.] | units.m
sun_and_stone[1].velocity = [500794.477878, 373808.365427, 0.] | (units.m/units.s)
sun_and_stone[1].mass = 0. | units.kg

print " particles:"
print sun_and_stone

converter=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
kepler=Kepler(converter)
kepler.initialize_code()
kepler.initialize_from_particles(sun_and_stone)
semi, ecc = kepler.get_elements()
apo = kepler.get_apastron()
per = kepler.get_periastron()
period = kepler.get_period()
kepler.stop()

print " semi-major axis ", semi.in_(units.AU)
print " eccentricity ", ecc
print " apocenter ", apo.in_(units.AU)
print " pericenter ", per.in_(units.AU)
print " period", period


