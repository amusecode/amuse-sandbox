import numpy 
numpy.random.seed(122222)

from amuse.units import units
from amuse.units import nbody_system
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from amuse.legacy.sse.interface import SSE
from amuse.legacy.mesa.interface import MESA


from amuse.datamodel import particle_attributes
from amuse.ic.plummer import new_plummer_sphere
from amuse.ic.salpeter import new_salpeter_mass_distribution
def demo3(N):

  salpeter_masses = new_salpeter_mass_distribution(N)
  total_mass = salpeter_masses.sum()

  convert_nbody = nbody_system.nbody_to_si(total_mass, 1.0 | units.parsec)

  print convert_nbody.to_nbody(100 | units.Myr)

  parts = new_plummer_sphere(N,convert_nbody)
  parts.radius=0. | units.RSun
  parts.mass = salpeter_masses
  parts.move_to_center()

  gravity=PhiGRAPE(convert_nbody,use_gl=True)

  eps=0.001 | units.parsec
  gravity.parameters.epsilon_squared = eps**2 

  gravity.particles.add_particles(parts)

  stellar_evolution = SSE()
  stellar_evolution.initialize_module_with_default_parameters() 

  stellar_evolution.particles.add_particles(parts)

  return gravity,stellar_evolution

def joined_evolve(gravity,stellar_evolution,end_time,dt):
  from_stellar_evolution_to_gravity = stellar_evolution.particles.new_channel_to(gravity.particles)
  time=gravity.model_time
  initial_mass=gravity.particles.mass.sum()
  while time < end_time:
    time += dt
    gravity.evolve_model(time)
    stellar_evolution.evolve_model(time)
    from_stellar_evolution_to_gravity.copy_attributes(["mass"])
    print 'mass fraction:',(gravity.particles.mass.sum()/initial_mass).number
    
    
if __name__=="__main__":
  gravity,stellar_evolution=demo3(100)
  gravity.start_viewer()
  joined_evolve(gravity,stellar_evolution, 100. | units.Myr , 0.25 | units.Myr)
