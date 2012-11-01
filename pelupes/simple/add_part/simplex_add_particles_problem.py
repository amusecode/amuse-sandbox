from amuse.community.simplex.interface import SimpleX
from amuse.units import units
from amuse.ext.molecular_cloud import ism_cube

from amuse.datamodel import Particle

from amuse.ext.evrard_test import body_centered_grid_unit_cube,sobol_unit_cube,regular_grid_unit_cube,uniform_random_unit_cube
import numpy
import pickle


N=500
Lstar=100. | units.LSun
L=10. | units.parsec
rho=1. | (units.amu/units.cm**3)
xion=0.0001 | units.none
u=(9. |units.kms)**2

Lsource=Lstar/ (20. | units.eV)

import logging
logging.basicConfig(level=logging.DEBUG)

print "Lsource:",Lsource.in_(units.s**-1)
print "timescale:",((rho*L**3/(1| units.amu) )/Lsource).in_(units.Myr)

rad = SimpleX()
rad.parameters.box_size=2.001*L    
rad.parameters.timestep=0.001 | units.Myr
#rad.parameters.number_of_border_sites=5000
print rad.parameters

particles=ism_cube(N,L/2.,rho,u,base_grid=uniform_random_unit_cube).result
particles.rho = rho
particles.flux = 0. | units.s**-1
particles.xion=xion

toadd=particles[0:100].copy()
particles=particles[100:].copy()


source=particles.add_particle(Particle())
source.x=0.|units.kpc
source.y=0.|units.kpc
source.z=0.|units.kpc
source.rho = rho 
source.flux = Lsource
source.xion = xion 
source.u = u

rad.particles.add_particles(particles)
rad.evolve_model(0.01 | units.yr)


bla = rad.particles.add_particles(toadd)
rad.recommit_particles()


print bla.xion

print len(rad.particles)


rad.stop()
