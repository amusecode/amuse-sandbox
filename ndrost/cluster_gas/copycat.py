from amuse.units import nbody_system
from amuse.units import units

class copycat(object):
  def __init__(self,baseclass, systems, converter, parameters=None, extra=dict()):
    self.baseclass=baseclass
    self.systems=systems
    self.converter=converter
    self.parameters=parameters
    self.extra=extra
        
  def get_gravity_at_point(self,radius,x,y,z):
    instance=self.baseclass(self.converter, **self.extra)
    for param,value in self.parameters:
      err=instance.parameters.__setattr__(param,value)
    for system in self.systems:
      instance.particles.add_particles(system.particles)
    ax,ay,az=instance.get_gravity_at_point(radius,x,y,z)
    instance.stop()
    return ax,ay,az

  def get_potential_at_point(self,radius,x,y,z):
    instance=self.baseclass(self.converter, **self.extra)
    for param,value in self.parameters:
      err=instance.parameters.__setattr__(param,value)
    for system in self.systems:
      instance.particles.add_particles(system.particles)
    phi=instance.get_potential_at_point(radius,x,y,z)
    instance.stop()
    return phi

class reinitializecopycat(object):
  def __init__(self,baseclass, systems, converter, parameters=None, extra=dict()):
    self.baseclass=baseclass
    self.systems=systems
    self.converter=converter
    self.parameters=parameters
    self.extra=extra
    self.instance=self.baseclass(self.converter, **self.extra)
        
  def get_gravity_at_point(self,radius,x,y,z):
    self.instance.initialize_code()
    for param,value in self.parameters:
      err=self.instance.parameters.__setattr__(param,value)
    for system in self.systems:
      self.instance.particles.add_particles(system.particles)
    ax,ay,az=self.instance.get_gravity_at_point(radius,x,y,z)
    self.instance.cleanup_code()
    return ax,ay,az

  def get_potential_at_point(self,radius,x,y,z):
    self.instance.initialize_code()
    for param,value in self.parameters:
      err=self.instance.parameters.__setattr__(param,value)
    for system in self.systems:
      self.instance.particles.add_particles(system.particles)
    phi=self.instance.get_potential_at_point(radius,x,y,z)
    self.instance.cleanup_code()
    return phi

class resetparticlescopycat(object):
  def __init__(self,baseclass, systems, converter, parameters=None, extra=dict()):
    self.baseclass=baseclass
    self.systems=systems
    self.converter=converter
    self.parameters=parameters
    self.extra=extra
    self.instance=self.baseclass(self.converter, **self.extra)
    self.instance.initialize_code()

  def get_gravity_at_point(self,radius,x,y,z):
    for param,value in self.parameters:
      err=self.instance.parameters.__setattr__(param,value)
    for system in self.systems:
      self.instance.particles.add_particles(system.particles)
    ax,ay,az=self.instance.get_gravity_at_point(radius,x,y,z)
    for system in self.systems:
        self.instance.particles.remove_particles(system.particles)
    return ax,ay,az

  def get_potential_at_point(self,radius,x,y,z):
    for param,value in self.parameters:
      err=self.instance.parameters.__setattr__(param,value)
    for system in self.systems:
      self.instance.particles.add_particles(system.particles)
    phi=self.instance.get_potential_at_point(radius,x,y,z)
    for system in self.systems:
        self.instance.particles.remove_particles(system.particles)
    return phi
