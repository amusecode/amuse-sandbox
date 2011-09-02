import numpy 
numpy.random.seed(122222)

from amuse.units import units
from amuse.units import nbody_system
from amuse.ext.plummer import MakePlummerModel
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from amuse.ext.salpeter import SalpeterIMF
from amuse.support.data import particle_attributes

def demo2(N):

  initial_mass_function = SalpeterIMF()
  total_mass, salpeter_masses = initial_mass_function.next_set(N)

  convert_nbody = nbody_system.nbody_to_si(total_mass, 1.0 | units.parsec)

  parts=MakePlummerModel(N,convert_nbody).result
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
