from amuse.units import units
from amuse.units import constants
from amuse.community.simplex.interface import SimpleX

from amuse.ext.evrard_test import uniform_random_unit_cube,uniform_unit_sphere

from matplotlib import pyplot

from amuse import datamodel
Tinit=100. | units.K
Tion=10000. | units.K
gamma1=1.666667-1
rhoinit=0.001 | (units.amu / units.cm**3)
uinit=constants.kB * Tinit / gamma1
uion=constants.kB * Tion / gamma1

def fake_internal_energy_from_xion(xion):
  return uinit+(uion-uinit)*xion 

def iliev_test_5( N=10000,
                  Ns=10,
                  L=15. | units.kpc ):

 
  x,y,z=uniform_random_unit_cube(N).make_xyz()
  p=datamodel.Particles(N)
  p.x=L*x+L
  p.y=L*y+L
  p.z=L*z+L
  p.vx= 0. | (units.km/units.s)
  p.vy= 0. | (units.km/units.s)
  p.vz= 0. | (units.km/units.s)
  p.u= uinit 
  p.rho = rhoinit
  p.flux=0. | (units.s**-1)
  p.xion=0. | units.none

  sources=datamodel.Particles(Ns)
  x,y,z=uniform_unit_sphere(Ns).make_xyz()

  sources.x=L*x*(1./N)**(1./3)+L
  sources.y=L*y*(1./N)**(1./3)+L
  sources.z=L*z*(1./N)**(1./3)+L
  sources.rho=0.001|(units.amu/units.cm**3)
  sources.flux=(5.e48/Ns) | (units.s**-1)
  sources.xion=1. | units.none
  
  return p,sources

interface=SimpleX(number_of_workers = 1,redirection='none')
interface.initialize_code()

L=15 | units.kpc
interface.parameters.box_size=2*L
interface.parameters.hilbert_order=1

p,s=iliev_test_5(N=1000,Ns=10,L=L)

p.add_particles(s)

interface.particles.add_particles(p)

interface.commit_particles()

interface.evolve_model( 0.05 | units.Myr)

p=interface.particles.copy()

cx=cy=cz=L
r=(((p.x-cx)**2+(p.y-cy)**2+(p.z-cz)**2) ** 0.5).value_in(units.kpc)
xion=p.xion.number
t=gamma1/constants.kB*fake_internal_energy_from_xion(xion)

t=t.value_in(units.K)

pyplot.figure(figsize=(8,6))
pyplot.semilogy(r,t,'r .')
pyplot.savefig("test.png")
