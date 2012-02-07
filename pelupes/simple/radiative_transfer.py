from amuse.community.simplex2_5.interface import SimpleX
from amuse.units import units
from amuse.ext.molecular_cloud import ism_cube

from amuse.datamodel import Particle

N=1000
L=13000. | units.parsec
rho=0.001 | (units.amu/units.cm**3)
xion=0.0001 | units.none
Lsource=5.0e48 | units.s**-1
u=(9. |units.kms)**2

rad = SimpleX()
        
particles=ism_cube(N,L/2,rho,u).result
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

rad.evolve_model(10. | units.Myr)

print "min Xion:", rad.particles.xion.amin()
print "average Xion:", rad.particles.xion.mean()
print "max Xion:", rad.particles.xion.amax()

rad.stop()
