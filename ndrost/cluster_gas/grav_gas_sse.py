"""
notes: contains ugly fix with time_offset to keep track of total time between 
restarts

version with feedback safety

"""
import numpy
    
from amuse.support.units import nbody_system
from amuse.support.units import units
    
from amuse import datamodel
from amuse.support.data.values import zero

from fast import FAST
from lmech import lmech
from copycat import copycat
from copycat import reinitializecopycat
from copycat import resetparticlescopycat
from amuse.ext.evrard_test import uniform_unit_sphere

numpy.random.seed(123456)

def sys_from_parts(base_class, parts=None, gasparts=None, parameters=None,converter=None, extra=dict()):
  sys=base_class(converter, **extra)
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
               conv,mgas,star_parts,gas_parts,dt_feedback, dt_fast,
               grav_parameters=(),gas_parameters=(),couple_parameters=(),
               feedback_efficiency=.01, feedback_radius=0.01 | units.parsec,
               total_feedback_energy=zero,evo_particles=None, star_particles_addition=None,
               start_time_offset=zero, feedback_safety=1.e-4*1.e51 | units.erg,
               feedback_dt=(10. | units.yr, 7500. | units.yr), feedback_period=20000. | units.yr,
               feedback_lasttime=None,
               
               grav_code_extra=dict(mode='gpu', redirection='none'),
               gas_code_extra=dict(number_of_workers=3,use_gl=False, redirection='none'),
               se_code_extra=dict(redirection='none'),
               grav_couple_code_extra=dict()
               ):
    
    self.codes=(gas_code,grav_code,se_code,grav_couple_code)
    self.parameters=(grav_parameters,gas_parameters,couple_parameters)
    self.conv=conv
    
    self.sph=sys_from_parts(gas_code, gasparts=gas_parts,
                         parameters=gas_parameters,
                         converter=conv,extra=gas_code_extra)                         
    self.grav=sys_from_parts(grav_code, parts=star_parts, 
                          parameters=grav_parameters,
                          converter=conv,extra=grav_code_extra)

    print gas_parameters
    print
    print self.sph.parameters
  
    self.sph_grav=reinitializecopycat(grav_couple_code, (self.sph,self.grav), conv,
                            parameters=couple_parameters, extra=grav_couple_code_extra)
    self.star_grav=reinitializecopycat(grav_couple_code, (self.sph,), conv,
                            parameters=couple_parameters, extra=grav_couple_code_extra)
  
    self.fast=FAST(verbose=True, timestep=dt_fast)
    self.fast.add_system(self.sph, (self.sph_grav,),False)
    self.fast.add_system(self.grav, (self.star_grav,),False)

    self.evo=se_code(**se_code_extra)
    self.evo.initialize_module_with_default_parameters() 
    if evo_particles is None:
      self.evo.particles.add_particles(star_parts)
    else:
      if len(evo_particles)!=len(star_parts):
        print "evo parts != star parts"
        raise Exception     
      self.evo.particles.add_particles(evo_particles)
      
    self.total_feedback_energy=total_feedback_energy
    self.dt_feedback=dt_feedback
    self.dt_fast=dt_fast
    self.mgas=mgas

    self.time=self.fast.model_time+start_time_offset
    self.time_offset=start_time_offset            
    if star_particles_addition is None:
      self.star_particles_addition=star_parts.empty_copy()
      self.star_particles_addition.Emech_last_feedback=0. | units.erg
    else:
      self.star_particles_addition=star_particles_addition

    self.feedback_efficiency=feedback_efficiency
    self.feedback_radius=feedback_radius
    
    self.feedback_safety=feedback_safety
    self.feedback_dt=feedback_dt    
    self.feedback_period=feedback_period
    if self.feedback_period < 1.0001*self.dt_feedback:
      print "feedback_period < dt_feedback, resetting to dt_feedback=", \
        self.dt_feedback.in_(units.yr)
      self.feedback_period=1.0001*self.dt_feedback
    if feedback_lasttime is None:
      self.feedback_lasttime=self.time-2*self.feedback_period
    else:
      self.feedback_lasttime=feedback_lasttime

  def evolve_model(self,tend):
    dt=self.dt_feedback
    time=self.time
    
    while time<tend-dt/2:
      time+=dt
      print self.feedback_lasttime.in_(units.yr) ,time.in_(units.yr) ,self.feedback_period.in_(units.yr) 
      if time > self.feedback_lasttime+self.feedback_period:
        self.sph.parameters.max_size_timestep=self.feedback_dt[1]
        print "long timestep"      
      else:
        self.sph.parameters.max_size_timestep=self.feedback_dt[0]
        print "short timestep"      

      self.fast.evolve_model(time-self.time_offset)
      print self.fast.model_time,',',self.sph.model_time,self.grav.model_time
      self.evo.evolve_model(self.fast.model_time)
      energy_added=self.mechanical_feedback()

      print energy_added.in_(units.erg), self.feedback_safety.in_(units.erg)
      print self.feedback_lasttime.in_(units.yr) ,time.in_(units.yr) 
      if energy_added > self.feedback_safety:
        self.feedback_lasttime=time
          
    self.time=self.fast.model_time+self.time_offset
            
  def mechanical_feedback(self):
    
    star_particles=self.star_particles.copy_to_memory()
    
    channel=self.particles.new_channel_to(star_particles)

    channel.copy_attributes(["x","y","z","vx","vy","vz"])
    channel.copy_attribute("mass","grav_mass")
    del channel

    channel=self.star_particles_addition.new_channel_to(star_particles)
    channel.copy_attribute("Emech_last_feedback")
    del channel

    new_sph=datamodel.Particles(0)
    
    star_particles.dmass=star_particles.grav_mass-star_particles.mass
    star_particles.u=(star_particles.Emech-star_particles.Emech_last_feedback)/star_particles.dmass    
    if numpy.any((star_particles.Emech-star_particles.Emech_last_feedback) < zero):
      print "feedback error"
      raise Exception
          
    losers=star_particles.select_array(lambda x: x > self.mgas, ["dmass"])
    while len(losers) > 0:            
      add=datamodel.Particles(len(losers))
      add.mass=self.mgas
      add.h_smooth=0. | units.parsec
      dx,dy,dz=uniform_unit_sphere(len(losers)).make_xyz()
      add.x=losers.x+self.feedback_radius*dx      
      add.y=losers.y+self.feedback_radius*dy
      add.z=losers.z+self.feedback_radius*dz
      add.vx=losers.vx      
      add.vy=losers.vy
      add.vz=losers.vz
      add.u=self.feedback_efficiency*losers.u
      losers.grav_mass-=self.mgas
      losers.Emech_last_feedback+=self.mgas*losers.u
      new_sph.add_particles(add)  
      losers=star_particles.select_array(lambda x,y: x-y > self.mgas, ["grav_mass","mass"])


    print "gas particles added:", len(new_sph)
    if len(new_sph) == 0:
      return zero     
    self.sph.gas_particles.add_particles(new_sph)
    feedback_energy_added=(new_sph.mass*new_sph.u).sum()
    self.total_feedback_energy=self.total_feedback_energy+feedback_energy_added

    channel = star_particles.new_channel_to(self.particles)
    channel.copy_attribute("grav_mass","mass")  
    del channel

    channel=star_particles.new_channel_to(self.star_particles_addition)
    channel.copy_attribute("Emech_last_feedback")
    del channel
    
    return feedback_energy_added

  def synchronize_model(self):
    return self.fast.synchronize_model()

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

  def dump_system_state(self,filename):    
    from amuse.support.io import write_set_to_file
    import cPickle
    write_set_to_file(self.grav.particles,filename+".grav","amuse",append_to_file=False)  
    write_set_to_file(self.gas_particles,filename+".gas","amuse",append_to_file=False)  
    write_set_to_file(self.star_particles,filename+".evo","amuse",append_to_file=False)  
    write_set_to_file(self.star_particles_addition,filename+".add","amuse",append_to_file=False)  
    f=open(filename+".data",'wb')
    print self.total_feedback_energy
    cPickle.dump((self.codes,
                  self.conv,
                  self.parameters,
                  self.mgas,
                  self.feedback_efficiency,
                  self.feedback_radius,
                  self.time,
                  self.dt_feedback,
                  self.dt_fast,
                  self.total_feedback_energy.in_(1.e51*units.erg),
                  self.feedback_safety,
                  self.feedback_dt,
                  self.feedback_period,
                  self.feedback_lasttime
                  ),f)
    f.close()

  @classmethod
  def load_system_state(cls,filename,new_gas_options=(),
               grav_code_extra=dict(mode='gpu', redirection='none'),
               gas_code_extra=dict(number_of_workers=3,use_gl=False, redirection='none'),
               se_code_extra=dict(redirection='none'),
               grav_couple_code_extra=dict()):    
    from amuse.support.io import read_set_from_file
    import cPickle
    star_parts=read_set_from_file(filename+".grav",'amuse')
    gas_parts=read_set_from_file(filename+".gas",'amuse')
    evo=read_set_from_file(filename+".evo",'amuse')
    add=read_set_from_file(filename+".add",'amuse')
    f=open(filename+".data",'r')
    data=cPickle.load(f)
    f.close()
    gas_code=data[0][0]
    grav_code=data[0][1]
    se_code=data[0][2]
    grav_couple_code=data[0][3]
    conv=data[1]
    gravp=data[2][0]
    gasp=data[2][1]
    cp=data[2][2]
    mgas=data[3]
    fe=data[4]
    fr=data[5]
    to=data[6]
    dt_feedback=data[7]
    dt_fast=data[8]
    tfe=data[9]
    fs=data[10]
    fdt=data[11]
    fp=data[12]
    flt=data[13]
    
    gasp=gasp+new_gas_options    
        
    return conv,cls(grav_code,gas_code,se_code,grav_couple_code, 
               conv,mgas,star_parts,gas_parts,dt_feedback, dt_fast,
               grav_parameters=gravp,gas_parameters=gasp,couple_parameters=cp,
               feedback_efficiency=fe, feedback_radius=fr,
               total_feedback_energy=tfe,evo_particles=evo, star_particles_addition=add,
               start_time_offset=to,feedback_safety=fs,
               feedback_dt=fdt, feedback_period=fp,
               feedback_lasttime=flt,
	       grav_code_extra=grav_code_extra,
               gas_code_extra=gas_code_extra,
               se_code_extra=se_code_extra,
               grav_couple_code_extra=grav_couple_code_extra)
