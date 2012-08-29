from amuse.community.fi.interface import Fi
from amuse.community.gadget2.interface import Gadget2
from amuse.datamodel import Particles

"""
class for radiative hydrodynamics

"""


class RadiativeHydro(object):
    
    def __init__(self, rad=None,hydro=None,src_particles=None):
        if rad is None:
          raise Exception("provide with rad class")
        if hydro is None:
          raise Exception("provide with hydro (sph) class")
        self.rad=rad()
        self.hydro=hydro()
        self._src_particles=src_particles

        self.hydro_parameters=self.hydro.parameters
        self.rad_parameters=self.rad.parameters
        self.rad_particles=self.rad.gas_particles
        self.src_particles=self.rad.src_particles
        self.gas_particles=self.hydro.gas_particles
        self.particles=self.hydro.particles
        self.dm_particles=self.hydro.dm_particles
        try:
          self.star_particles=self.hydro.star_particles
        except:
          self.star_particles=None
        self.timestep=None

    @property
    def model_time(self):
        return self.hydro.model_time

    def evolve_model(self,tend):
  
        if self._src_particles:
          self._src_particles.synchronize_to(self.rad.src_particles)
          self._src_particles.new_channel_to(self.rad.sr_particles).copy_attributes(["x","y","z","luminosity"])
  
        self.hydro.gas_particles.synchronize_to(self.rad.gas_particles)
    
        update_rad_from_hydro_channel = self.hydro.gas_particles.new_channel_to(self.rad.gas_particles)
        update_hydro_from_rad_channel = self.rad.gas_particles.new_channel_to(self.hydro.gas_particles)

        rad_attributes_to_update=["x","y","z","rho","h_smooth"]
        hydro_attributes_to_update=["u"]
        try:
          momentum_kicks=self.rad.parameters.momentum_kicks_flag
          if momentum_kicks:
            rad_attributes_to_update.extend(["vx","vy","vz"])
            hydro_attributes_to_update.extend(["vx","vy","vz"])
        except:
          momentum_kicks=False
    
        t=self.model_time
        timestep=self.timestep
        if timestep is None:
          timestep=tend-t
    
        while t < (tend-timestep/2):
          self.hydro.evolve_model(t+timestep/2.)
          update_rad_from_hydro_channel.copy_attributes(rad_attributes_to_update)
          self.rad.evolve_model(t+timestep)
          update_hydro_from_rad_channel.copy_attributes(hydro_attributes_to_update)
          self.hydro.evolve_model(t+timestep)
          t=self.model_time
    
    def radhydro_particles_copy(self):
          parts=self.rad.gas_particles.copy()
          channel = self.gas_particles.new_channel_to(parts)
          channel.copy_attributes(["mass","x","y","z","vx","vy","vz","rho","u"])  	  
          return parts

