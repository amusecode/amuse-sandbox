import numpy 
numpy.random.seed(122222)

from amuse.units import nbody_system
from amuse.ext.plummer import MakePlummerModel
from amuse.legacy.phiGRAPE.interface import PhiGRAPE


def demo1(N):
  parts=MakePlummerModel(N).result
  parts.radius=0. | nbody_system.length

  #interface=PhiGRAPE(use_gl=True)
  interface=PhiGRAPE()	
  eps=0.001 | nbody_system.length
  interface.parameters.epsilon_squared = eps**2 

  interface.particles.add_particles(parts)

  return interface

if __name__=="__main__":
  interface=demo1(100)
  interface.start_viewer()
  interface.evolve_model(100. | nbody_system.time)
