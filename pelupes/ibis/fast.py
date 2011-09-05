# issues:
# - for now, units in si 
# - a common coordinate system is used for all systems
# - sync of systems should be checked
# - timestepping: adaptive dt?

from amuse.units import units
import threading

from amuse.support import data
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
#  parts.copy_values_of_state_attributes_to(system.particles)
  channel = parts.new_channel_to(system.particles)
  channel.copy_attributes(["vx","vy","vz"])
  
class FAST(object):
  def __init__(self,verbose=False):
    """
    verbose indicates whether to output some run info
    """  
    self.systems=set()
    self.partners=dict()
    self.time_offsets=dict()
    self.time=0. | units.s
    self.do_sync=dict()
    self.verbose=verbose
  
  def add_system(self, interface,  partners=set(),do_sync=True):
    """
    add a system to bridge integrator  
    """
    if hasattr(interface,"model_time"):
      self.time_offsets[interface]=(self.time-interface.model_time)
    else:
      self.time_offsets[interface]=0.        
    self.systems.add(interface)
    for p in partners:
      if not hasattr(interface,"get_gravity_at_point"):
          return -1
    self.partners[interface]=partners
    self.do_sync[interface]=do_sync  
    return 0
    
  def evolve_model(self,tend,timestep=None):
    """
    evolve combined system to tend, timestep fixes timestep
    """
    if timestep is None:
      timestep=tend-self.time
    while self.time < tend:    
      dt=min(timestep,tend-self.time)
      dt=timestep
      self.kick_systems(dt/2)   
      self.drift_systems(self.time+dt)
      self.kick_systems(dt/2)
      self.time=self.time+dt
    return 0    
  
  def synchronize_model(self):
    """ 
    explicitly synchronize all components
    """
    for x in self.systems:
      if hasattr(x,"synchronize_model"):
        if(self.verbose): print x.__class__.__name__,"is synchronizing",
        x.synchronize_model()    
        if(self.verbose): print ".. done"
                          
  def get_potential_at_point(self,radius,x,y,z):
    err=0
    pot=0.*radius
    for x in self.systems:
      _pot,err=x.get_potential_at_point(radius,x,y,z)
      if err != 0: 
        break
      pot=pot+_pot
    return pot,err
      
  def get_gravity_at_point(self,radius,x,y,z):
    err=0
    ax=0.*radius
    ay=0.*radius
    az=0.*radius
    for x in self.systems:
      _ax,_ay,_az,err=x.get_gravity_at_point(radius,x,y,z)
      if err != 0: 
        break
      ax=ax+_ax
      ay=ay+_ay
      az=az+_az
    return ax,ay,az,err
    
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
  def thermal_energy(self):  
    Eth=0. | units.kg* (units.m/units.s)**2
    for x in self.systems:
      if hasattr(x,'thermal_energy'):
        Eth+=x.thermal_energy
    return Eth

  @property
  def model_time(self):  
    return self.time


        
  @property
  def particles(self):
    arr=[]
    for x in self.systems:
      if hasattr(x,"particles"):
        arr.append(x.particles)
    return data.ParticlesSuperset(arr)          

  @property
  def gas_particles(self):
    arr=[]
    for x in self.systems:
      if hasattr(x,"gas_particles"):
        arr.append(x.gas_particles)
    return data.ParticlesSuperset(arr)          

# 'private' functions

  def drift_systems(self,tend):
    threads=[]
    for x in self.systems:
      if hasattr(x,"evolve_model"):
        offset=self.time_offsets[x]
        if(self.verbose):
          print "evolving", x.__class__.__name__,
        threads.append(threading.Thread(target=x.evolve_model, args=(tend-offset,)) )
    for x in threads:
      x.start()
    for x in threads:
      x.join()
    if(self.verbose): 
      print ".. done"
    return 0

  def kick_systems(self,dt):
    for x in self.systems:
      if self.do_sync[x]:
        if hasattr(x,"synchronize_model"):
          if(self.verbose): print x.__class__.__name__,"is synchronizing",
          x.synchronize_model()    
          if(self.verbose):  print ".. done"
    for x in self.systems:
      if hasattr(x,"particles"):
        for y in self.partners[x]:
          if x is not y:
            if(self.verbose):  print x.__class__.__name__,"receives kick from",y.__class__.__name__,
            kick_system(x,y.get_gravity_at_point,dt)
            if(self.verbose):  print ".. done"
    return 0
