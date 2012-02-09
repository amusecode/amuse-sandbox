from amuse.community.simplex.interface import SimpleX
from amuse.units import units
from amuse.ext.molecular_cloud import ism_cube

from amuse.datamodel import Particle

N=1000
Lstar=100. | units.LSun
L=10. | units.parsec
rho=1. | (units.amu/units.cm**3)
xion=0.0001 | units.none
u=(9. |units.kms)**2

Lsource=Lstar/ (20. | units.eV)

print Lsource.in_(units.s**-1)
print ((rho*L**3/(1| units.amu) )/Lsource).in_(units.Myr)

rad = SimpleX()
rad.parameters.box_size=1.001*L    
rad.parameters.timestep=0.001 | units.Myr

particles=ism_cube(N,L/2.,rho,u).result
particles.rho = rho
particles.flux = 0. | units.s**-1
particles.xion=xion

source=particles.add_particle(Particle())
source.x=0.|units.kpc
source.y=0.|units.kpc
source.z=0.|units.kpc
source.rho = rho 
source.flux = Lsource
source.xion = xion 
source.u = u

rad.particles.add_particles(particles)

rad.evolve_model(0.1 | units.Myr)

print "min Xion:", rad.particles.xion.amin()
print "average Xion:", rad.particles.xion.mean()
print "max Xion:", rad.particles.xion.amax()

rad.stop()
