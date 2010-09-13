import numpy
    
from amuse.support.units import nbody_system
from amuse.support.units import units
    
from amuse.support.data import core
from amuse.support.data.values import zero

from fast import FAST
from lmech import lmech
from lmech import e_supernova
from copycat import copycat
from amuse.ext.evrard_test import uniform_unit_sphere

numpy.random.seed(123456)

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

class grav_gas_sse(object):
  def __init__(self,grav_code,gas_code,se_code,grav_couple_code, 
               conv,mgas,star_parts,gas_parts,eps,dt_feedback, dt_fast,
               grav_parameters=(),gas_parameters=(),fast_parameters=(),
               feedback_efficiency=.01, feedback_radius=0.01 | units.parsec):
    
    self.sph=sys_from_parts(gas_code, gasparts=gas_parts,
                         parameters=gas_parameters,
                         converter=conv,extra=dict(use_gl=False))
    self.grav=sys_from_parts(grav_code, parts=star_parts, 
                          parameters=grav_parameters,
                          converter=conv,extra=dict(use_gl=False))
 
    self.sph_grav=copycat(grav_couple_code, (self.sph,self.grav), conv,eps2=eps**2)
    self.star_grav=copycat(grav_couple_code, (self.sph,), conv,eps2=eps**2)
  
    self.fast=FAST(verbose=True, timestep=dt_fast)
    self.fast.add_system(self.sph, (self.sph_grav,),False)
    self.fast.add_system(self.grav, (self.star_grav,),False)

    self.evo=se_code()
    self.evo.initialize_module_with_default_parameters() 
    self.evo.particles.add_particles(star_parts)

    self.total_feedback_energy=zero
    self.dt_feedback=dt_feedback
    self.dt_fast=dt_fast
    self.mgas=mgas

    self.time=self.fast.model_time        
    self.star_particles_addition=star_parts.empty_copy()
    self.star_particles_addition.time_last_feedback=self.time

    self.feedback_efficiency=feedback_efficiency
    self.feedback_radius=feedback_radius

  def evolve_model(self,tend):
    dt=self.dt_feedback
    time=self.time
    
    while (time<tend-dt/2):
      time+=dt
      self.fast.evolve_model(time)
      print self.fast.model_time,',',self.sph.model_time,self.grav.model_time
      self.evo.evolve_model(self.fast.model_time)
#      self.mechanical_feedback(self.fast.model_time)
    self.time=self.fast.model_time        
#    print self.time
#    raise Exception
            
  def mechanical_feedback(self,time):
    
    star_particles=self.star_particles.copy_to_memory()
    
    channel=self.particles.new_channel_to(star_particles)
    channel.copy_attribute("mass","grav_mass")
    del channel

    channel=self.star_particles_addition.new_channel_to(star_particles)
    channel.copy_attribute("time_last_feedback")
    del channel
    
    new_sph=core.Particles(0)
    
    star_particles.dtime=time-star_particles.time_last_feedback
    star_particles.dmass=star_particles.grav_mass-star_particles.mass
    losers=star_particles.select_array(lambda x,y: x-y > self.mgas, ["grav_mass","mass"])
    losers.time_last_feedback=time
    print losers.mass.value_in(units.MSun)
    while len(losers) > 0:
      add=core.Particles(len(losers))
      add.mass=self.mgas
      add.h_smooth=0. | units.parsec
      dx,dy,dz=uniform_unit_sphere(len(losers)).make_xyz()
      add.x=losers.x+self.feedback_radius*dx      
      add.y=losers.y+self.feedback_radius*dy
      add.z=losers.z+self.feedback_radius*dz
      add.vx=losers.vx      
      add.vy=losers.vy
      add.vz=losers.vz
      add.u=self.feedback_efficiency*(lmech(losers)* \
              losers.dtime/losers.dmass+e_supernova(losers)/losers.dmass)
      losers.grav_mass-=self.mgas
      new_sph.add_particles(add)  
      losers=star_particles.select_array(lambda x,y: x-y > self.mgas, ["grav_mass","mass"])


    print "gas particles added:", len(new_sph)
    if len(new_sph) == 0:
      return      
    self.sph.gas_particles.add_particles(new_sph)
    self.total_feedback_energy=self.total_feedback_energy+(new_sph.mass*new_sph.u).sum()
    print new_sph.u**0.5
    
    channel = star_particles.new_channel_to(self.particles)
    channel.copy_attribute("grav_mass","mass")  
    del channel

    channel=star_particles.new_channel_to(self.star_particles_addition)
    channel.copy_attribute("time_last_feedback")
    del channel

  @property
  def kinetic_energy(self):
    return self.fast.kinetic_energy
 
  @property     
  def potential_energy(self):
    return self.fast.potential_energy
  
  @property
  def thermal_energy(self):
    return self.fast.thermal_energy
  
  @property
  def feedback_energy(self):
    return self.total_feedback_energy
  
  @property
  def model_time(self):
    return self.time
  
  @property
  def particles(self):
    return self.fast.particles
      
  @property
  def gas_particles(self):
    return self.fast.gas_particles
    
  @property
  def star_particles(self):
    return self.evo.particles

    
  
