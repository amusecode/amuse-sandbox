import cPickle

from amuse.community.huayno.interface import Huayno
from amuse.community.sse.interface import SSE

from amuse.datamodel import Particles
from amuse.units import units

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
      if not kwargs.has_key("begin_time"):
        kwargs["begin_time"]=self.grav.model_time
      self.grav.parameters.begin_time=kwargs['begin_time']
      self.time=kwargs['begin_time']
      self.timestep=kwargs['timestep']
      self.particles_accessed=False
    @property
    def particles(self):
      return self._particles
    def evolve_model(self,tend):
        
      self._particles.synchronize_to(self.grav.particles)
      channel=self._particles.new_channel_to(self.grav.particles)
      channel.copy_attributes(["mass","radius","x","y","z","vx","vy","vz"])

      if hasattr(self._particles,"is_star"):
        stars=self._particles.select_array(lambda x:x, ["is_star"]).copy()
      else:
        stars=self._particles.copy()

      remove_set=self.se.particles.difference(stars)
      if len(remove_set)>0:        
        self.se.particles.remove_particles(remove_set)

      add_set=stars.difference(self.se.particles).copy()

      if len(add_set)>0: 
        if hasattr(add_set,"zams_mass"):
          add_set.mass=add_set.zams_mass  
        self.se.particles.add_particles(add_set)
  
      channel=self.se.particles.new_channel_to(self.grav.particles)
      while self.time<(tend-self.timestep/2):
        self.grav.evolve_model(self.time+self.timestep/2)
        self.se.evolve_model(self.time+self.timestep)
        channel.copy_attribute("mass")
        self.grav.evolve_model(self.time+self.timestep)
        self.time+=self.timestep
      
      print self._particles.mass[0].in_(units.MSun)
      
      channel=self.grav.particles.new_channel_to(self._particles)
      channel.copy_attributes(["mass","x","y","z","vx","vy","vz"])        

      channel=self.se.particles.new_channel_to(self._particles)
      attributes=set(self.se.particles.get_attribute_names_defined_in_store()) -\
                 set(["mass","x","y","z","vx","vy","vz"])
      channel.copy_attributes(attributes)
    @property  
    def model_time(self):  
      return self.time

if __name__=="__main__":

  from massloss import grav_with_se
  from remote_se import forwarding_class_client
  
  from amuse.ic.plummer import new_plummer_model
  
  from amuse.units import units,nbody_system
  
  print "START grav_w_se"
  
  N=10
  
  conv=nbody_system.nbody_to_si(N*22. | units.MSun,0.1| units.parsec)
  
  timestep=conv.to_si(0.25 | nbody_system.time)

  print timestep.in_(units.Myr)

#  code=grav_with_se(converter=conv)
  code=forwarding_class_client(grav_with_se,code_kwarg=dict(converter=conv,timestep=timestep))

  print code.parameters["SSE"]
  print code.parameters["Huayno"]
  
  p=new_plummer_model(N,conv)
  
  print p.mass.in_(units.MSun)
  
  code.particles.add_particles(p)
  
  code.evolve_model(10. | units.Myr)
  
  print code.particles.stellar_type
  print code.model_time
  
  
