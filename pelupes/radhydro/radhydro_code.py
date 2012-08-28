from amuse.community.fi.interface import Fi
from amuse.community.gadget2.interface import Gadget2
from amuse.datamodel import Particles

def radhydro_class(SPH):
  class radhydro(SPH):
    
    def __init__(self, *args, **kargs):
      if not kargs.has_key('rad'):
        raise Exception("provide with rad class factory")
      self._rad=kargs['rad']()
      if kargs.has_key('src_particles'):
        self._src_particles=src_particles
      else:
        self._src_particles=None
      SPH.__init__(self, *args, **kargs)
      self.hydro=self
      self.rad=self._rad
      self.hydro_parameters=self.parameters
      self.rad_parameters=self._rad.parameters
      self.rad_particles=self._rad.gas_particles
      self.src_particles=self._rad.src_particles

    def evolve_model(self,tend):
  
      if self._src_particles:
        self._src_particles.synchronize_to(self._rad.src_particles)
        self._src_particles.new_channel_to(self._rad.sr_particles).copy_attributes(["x","y","z","luminosity"])

      self.gas_particles.synchronize_to(self._rad.gas_particles)
  
      update_rad_from_hydro_channel = self.gas_particles.new_channel_to(self._rad.gas_particles)
      update_hydro_from_rad_channel = self.rad.gas_particles.new_channel_to(self.gas_particles)
  
      t=self.model_time
      timestep=self.timestep
      if timestep is None:
        timestep=tend-t
  
      while t < (tend-timestep/2):
        self.overridden().evolve_model(t+timestep/2.)
        update_rad_from_hydro_channel.copy_attributes(["x","y","z","rho","h_smooth"])
        self._rad.evolve_model(t+dt)
        update_hydro_from_rad_channel.copy_attributes(["u"])
        self.overridden().evolve_model(t+dt)
        t=self.model_time
  
    def radhydro_particles_copy(self):
        parts=self._rad.gas_particles.copy()
        channel = self.gas_particles.new_channel_to(parts)
        channel.copy_attributes(["mass","x","y","z","vx","vy","vz","rho","u"])  	  
        return parts

  return radhydro

fi_radhydro=radhydro_class(Fi)
gadget_radhydro=radhydro_class(Gadget2)

