# issues die nog nader onderzocht moeten worden:
# - for now, units in si 
# - init van posities, snelheden (relatief, absoluut?), precisie?
# - sync van systemen
# - timestepping: adaptive dt?

from amuse.support.data.values import zero
from amuse.support.units import units
from amuse.support.data import core

def potential_energy(system, get_potential):
  parts=system.particles.copy()
  pot=get_potential(parts.radius,parts.x,parts.y,parts.z)
  return (pot*parts.mass).sum()/2 

def kick_system(system, get_gravity, dt):
  parts=system.particles.copy()
  ax,ay,az=get_gravity(parts.radius,parts.x,parts.y,parts.z)
  parts.vx=parts.vx+dt*ax
  parts.vy=parts.vy+dt*ay
  parts.vz=parts.vz+dt*az
  parts.copy_values_of_all_attributes_to(system.particles)
  
class bridge(object):
  def __init__(self, do_sync=True, timestep=None):
    self.systems=set()
    self.partners=dict()
    self.time_offsets=dict()
    self.time=0. | units.s
    self.do_sync=do_sync
    self.timestep=timestep
  
  def add_system(self, interface, partners=set()):
    if hasattr(interface,"model_time"):
      self.time_offsets[interface]=(self.time-interface.model_time)
    else:
      self.time_offsets[interface]=zero
    self.systems.add(interface)
    for p in partners:
      if not hasattr(interface,"get_gravity_at_point"):
          return -1
    self.partners[interface]=partners  
    return 0
    
  def evolve_model(self,tend,timestep=None):
    if self.timestep is None:
      local_timestep=tend-self.time
    else:
      local_timestep=self.timestep
    while self.time < tend-local_timestep/2:    
      dt=min(local_timestep,tend-self.time)
      self.kick_systems(dt/2)   
      self.drift_systems(self.time+dt)
      self.kick_systems(dt/2)
      self.time=self.time+dt
    return 0    

  def drift_systems(self,tend):
    for x in self.systems:
      if hasattr(x,"evolve_model"):
        offset=self.time_offsets[x]
        print "evolving", x.__class__.__name__,
        x.evolve_model( tend-offset)
        print ".. done"
    return 0
    
  def kick_systems(self,dt):
    if(self.do_sync):
      for x in self.systems:
        if hasattr(x,"synchronize_model"):
          print x.__class__.__name__,"is synchronizing",
          x.synchronize_model()    
          print ".. done"
    for x in self.systems:
      if hasattr(x,"particles"):
        for y in self.partners[x]:
          if x is not y:
            print x.__class__.__name__,"receives kick from",y.__class__.__name__,
            kick_system(x,y.get_gravity_at_point,dt)
            print ".. done"
    return 0
                    
  def get_potential_at_point(self,radius,x,y,z):
    pot=zero
    for s in self.systems:
      _pot=s.get_potential_at_point(radius,x,y,z)
      pot=pot+_pot
    return pot
      
  def get_gravity_at_point(self,radius,x,y,z):
    ax=zero
    ay=zero
    az=zero
    for s in self.systems:
      _ax,_ay,_az=s.get_gravity_at_point(radius,x,y,z)
      ax=ax+_ax
      ay=ay+_ay
      az=az+_az
    return ax,ay,az
    
  @property
  def potential_energy(self):
    Ep=0. | units.kg* (units.m/units.s)**2
    for x in self.systems:
      Ep+=x.potential_energy
      if hasattr(x,"particles"):
        for y in self.partners[x]:
          Ep+=potential_energy(x,y.get_potential_at_point)
    return Ep
  
  @property
  def kinetic_energy(self):  
    Ek=0. | units.kg* (units.m/units.s)**2
    for x in self.systems:
      Ek+=x.kinetic_energy
    return Ek
        
  @property
  def particles(self):
    arr=[]
    for x in self.systems:
      if hasattr(x,"particles"):
        arr.append(x.particles)
    return core.ParticlesSuperset(arr)          
  
