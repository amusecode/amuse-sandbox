import numpy 
numpy.random.seed(122222)

from amuse.units import units
from amuse.units import nbody_system
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from amuse.datamodel import particle_attributes
from amuse.ic.plummer import new_plummer_sphere
from amuse.ic.salpeter import new_salpeter_mass_distribution
def demo2(N):

  salpeter_masses = new_salpeter_mass_distribution(N)
  total_mass = salpeter_masses.sum()

  convert_nbody = nbody_system.nbody_to_si(total_mass, 1.0 | units.parsec)

  parts = new_plummer_sphere(N,convert_nbody)
  parts.radius=0. | nbody_system.length
  parts.mass = salpeter_masses
  parts.move_to_center()

  interface=PhiGRAPE(convert_nbody,use_gl=True)

  eps=0.001 | units.parsec
  interface.parameters.epsilon_squared = eps**2 

  interface.particles.add_particles(parts)

  return interface

if __name__=="__main__":
  interface=demo2(100)
  interface.start_viewer()
  total_energy_0 = interface.particles.kinetic_energy() + interface.particles.potential_energy()   
  interface.evolve_model(100. | nbody_system.time)
  total_energy = interface.kinetic_energy + interface.potential_energy   
  print 'energy error:', ((total_energy-total_energy_0)/total_energy_0).number
