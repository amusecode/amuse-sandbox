from amuse.community.simplex.interface import SimpleX
from amuse.units import units
from amuse.ext.molecular_cloud import ism_cube

from amuse.datamodel import Particle

from amuse.ext.evrard_test import body_centered_grid_unit_cube,sobol_unit_cube,regular_grid_unit_cube,uniform_random_unit_cube
import numpy
import pickle

use_state_from_file = True
if use_state_from_file:
    with open('state.random','r') as stream:
        state = pickle.load(stream)
    numpy.random.set_state(state)
else:
    state = numpy.random.get_state()
    with open('state.random','w') as stream:
        pickle.dump(state, stream)


N=1000
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
rad.parameters.box_size=1.001*L    
rad.parameters.timestep=0.001 | units.Myr
#rad.parameters.hilbert_order=1
#rad.parameters.recombination_radiation_flag=0
#rad.parameters.number_of_freq_bins=5
#rad.parameters.thermal_evolution_flag=1
#rad.parameters.blackbody_spectrum_flag=1
#rad.parameters.metal_cooling_flag=0
print rad.parameters

#particles=ism_cube(N,L/2.,rho,u,base_grid=sobol_unit_cube).result
particles=ism_cube(N,L/2.,rho,u,base_grid=uniform_random_unit_cube).result
particles.rho = rho
particles.flux = 0. | units.s**-1
particles.xion=xion

toadd=particles[0:5].copy()
particles=particles[5:].copy()

print toadd

source=particles.add_particle(Particle())
source.x=0.|units.kpc
source.y=0.|units.kpc
source.z=0.|units.kpc
source.rho = rho 
source.flux = Lsource
source.xion = xion 
source.u = u

rad.particles.add_particles(particles)
rad.evolve_model(0.05 | units.Myr)

#print len(rad.particles)
#print rad.model_time.in_(units.Myr)
#print "min Xion:", rad.particles.xion.min()
#print "average Xion:", rad.particles.xion.mean()
#print "max Xion:", rad.particles.xion.max()
#print "max flux:", rad.particles.flux.max()

bla = rad.particles.add_particles(toadd)
print bla.index_in_code
rad.recommit_particles()
#print len(rad.particles)
#print rad.model_time.in_(units.Myr)
#print "min Xion:", rad.particles.xion.min()
#print "average Xion:", rad.particles.xion.mean()
#print "max Xion:", rad.particles.xion.max()
#print "max flux:", rad.particles.flux.max()

#print toadd.get_intersecting_subset_in(rad.particles).xion
print bla.xion
print bla.index_in_code
print len(rad.particles)


rad.stop()
