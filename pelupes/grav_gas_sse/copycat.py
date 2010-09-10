from amuse.support.units import nbody_system
from amuse.support.units import units


class copycat(object):
  def __init__(self,baseclass, systems, converter,eps2=0. | nbody_system.length**2):
    self.baseclass=baseclass
    self.systems=systems
    self.converter=converter
    self.epsilon_squared=eps2
        
  def get_gravity_at_point(self,radius,x,y,z):
    instance=self.baseclass(self.converter)
    instance.parameters.epsilon_squared = self.epsilon_squared
    for system in self.systems:
      instance.particles.add_particles(system.particles)
    ax,ay,az=instance.get_gravity_at_point(radius,x,y,z)
    instance.stop()
    return ax,ay,az

  def get_potential_at_point(self,radius,x,y,z):
    instance=self.baseclass(self.converter)
    instance.parameters.epsilon_squared = self.epsilon_squared
    for system in self.systems:
      instance.particles.add_particles(system.particles)
    phi=instance.get_potential_at_point(radius,x,y,z)
    instance.stop()
    return phi

def sys_from_parts(base_class, parts=None, gasparts=None, parameters=None,converter=None, extra=dict()):
  sys=base_class(convert_nbody=converter, **extra)
  sys.initialize_code()
  for param,value in parameters:
    err=sys.parameters.__setattr__(param,value)
  sys.commit_parameters()
  if parts is not None:
    sys.particles.add_particles(parts)
  if gasparts is not None:
    sys.gas_particles.add_particles(gasparts)
  if hasattr(sys,"start_viewer"):
    sys.start_viewer()
  return sys
