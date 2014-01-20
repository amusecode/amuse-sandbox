import cPickle

from amuse.community.huayno.interface import Huayno
from amuse.community.sse.interface import SSE

from amuse.datamodel import Particles

class combined_parameters(dict):
    def __init__(self, codes):
      for code in codes:
        self[code.__class__.__name__]=code.parameters
    def copy(self):
      copy=combined_parameters(())
      for code,p in self.items():
        copy[code]=p.copy()
      return copy
    def reset_from_memento(self,mementos):
      for code,p in mementos.items():
        self[code].reset_from_memento(p)

class grav_with_se(object):
    def __init__(self,*args,**kwargs):
      self.se=SSE(channel_type="sockets")
      self.grav=Huayno(kwargs["converter"],channel_type="sockets")
      self.parameters=combined_parameters((self.se,self.grav))
      self._particles=Particles()
      self.time=self.grav.model_time
      self.timestep=kwargs['timestep']
    @property
    def particles(self):
      return self._particles
    def evolve_model(self,tend):
      self._particles.synchronize_to(self.grav.particles)
      channel=self._particles.new_channel_to(self.grav.particles)
      channel.copy_all_attributes()

      self._particles.synchronize_to(self.se.particles)
      channel=self._particles.new_channel_to(self.se.particles)
      channel.copy_all_attributes()
      
      channel=self.se.particles.new_channel_to(self.grav.particles)
      while self.time<(tend-self.timestep/2):
        self.grav.evolve_model(self.time+self.timestep/2)
        channel.copy_attribute("mass")
        self.se.evolve_model(self.time+self.timestep)
        self.grav.evolve_model(self.time+self.timestep)
        self.time+=self.timestep
      
      channel=self.grav.particles.new_channel_to(self._particles)
      channel.copy_attributes(["x","y","z","vx","vy","vz"])        

      channel=self.se.particles.new_channel_to(self._particles)
      channel.copy_all_attributes()        
    @property  
    def model_time(self):  
      return self.time

#from massloss import grav_with_se

if __name__=="__main__":

  from massloss import grav_with_se
  from remote_se import forwarding_class_client
  
  from amuse.ic.plummer import new_plummer_model
  
  from amuse.units import units,nbody_system
  
  N=10
  
  conv=nbody_system.nbody_to_si(N*22. | units.MSun,0.1| units.parsec)
  
  timestep=conv.to_si(0.25 | nbody_system.time)

  print timestep.in_(units.Myr)

#  code=grav_with_se(converter=conv)
  code=forwarding_class_client(grav_with_se,code_kwarg=dict(converter=conv,timestep=timestep))
  
  p=new_plummer_model(N,conv)
  
  print p.mass.in_(units.MSun)
  
  code.particles.add_particles(p)
  
  code.evolve_model(10. | units.Myr)
  
  print code.particles.stellar_type
  print code.model_time
  
  
