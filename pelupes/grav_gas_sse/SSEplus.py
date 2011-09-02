import numpy

from amuse.community.sse.interface import SSE
from amuse.support.data import core

from amuse.units import units
from amuse.units.quantities import zero

from lmech import lmech

def e_supernova(stellar_type,prev_stellar_type):
  i=numpy.where( (stellar_type>=13 | units.stellar_type) &
                 (stellar_type<=15 | units.stellar_type) &
                 ((prev_stellar_type<13 | units.stellar_type) |
                  (prev_stellar_type>15 | units.stellar_type)) )[0] 
  n=len(stellar_type)
  e=numpy.array([0.]*n)
  e[i]=1.e51
  if(n == 1 ):
    return (units.erg).new_quantity(e[0])
  else:
    return (units.erg).new_quantity(e)

    
class SSEplus(SSE):
  def __init__(self,**options):
    self.model_time=zero
    SSE.__init__(self,**options)

  def evolve_model(self,tend):
    if not hasattr(self.particles,'Emech'):
      self.particles.Lmech=lmech(self.particles)
      self.particles.Emech=(0.| units.Myr)*self.particles.Lmech

    stellar_type=self.particles.stellar_type.copy()
    prev_lm=self.particles.Lmech.copy()
    
    ret=SSE.evolve_model(self,tend)
    if tend>self.model_time:
      dt=tend-self.model_time    
      self.model_time=tend
      lm=lmech(self.particles)
      self.particles.Lmech=lm.copy()
      self.particles.Emech=self.particles.Emech+dt*(prev_lm+lm)/2. + \
                          e_supernova(self.particles.stellar_type,stellar_type)    
    return ret        


if __name__=="__main__":
  evo=SSEplus()
  evo.initialize_module_with_default_parameters() 
  p=core.Particles(4)
  p.mass=units.MSun([15,25,35,50])
  evo.particles.add_particles(p)

  t=zero
  while t< 1 | units.Myr:
    t+=.125| units.Myr
    evo.evolve_model(t)
    
  print evo.particles.Emech.in_(1.e51*units.erg)
  print evo.model_time
