#import cProfile

import os
import sys
import numpy


try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    
from amuse.units import nbody_system
from amuse.units import units
    
from amuse.community.fi.interface import Fi
from amuse.community.bhtree.interface import BHTree

from amuse.ext.evrard_test2 import MakeEvrardTest
from amuse.ext.evrard_test2 import regular_grid_unit_cube
from amuse.ext.evrard_test2 import body_centered_grid_unit_cube
from amuse.test.amusetest import get_path_to_results

from fast import FAST
#from amuse.ext.derived_grav_systems import copycat
 
import logging
#logging.basicConfig(level=logging.DEBUG)


class copycat(object):
  def __init__(self,baseclass, system,converter,eps2=0. | nbody_system.length**2):
    self.baseclass=baseclass
    self.system=system
    self.converter=converter
    self.epsilon_squared=eps2
        
  def get_gravity_at_point(self,radius,x,y,z):
    instance=self.baseclass(self.converter,channel_type='sockets')
    instance.parameters.epsilon_squared = self.epsilon_squared
    instance.particles.add_particles(self.system.particles)
    ax,ay,az=instance.get_gravity_at_point(radius,x,y,z)
    instance.stop()
    return ax,ay,az

  def get_potential_at_point(self,radius,x,y,z):
    instance=self.baseclass(self.converter,channel_type='sockets')
    instance.parameters.epsilon_squared = self.epsilon_squared
    instance.particles.add_particles(self.system.particles)
    phi=instance.get_potential_at_point(radius,x,y,z)
    instance.stop()
    return phi


class copycat2(object):
  def __init__(self,baseclass, system,converter,eps2=0. | nbody_system.length**2):
    self.baseclass=baseclass
    self.system=system
    self.converter=converter
    self.epsilon_squared=eps2
    self.instance=self.baseclass(self.converter,channel_type='sockets')#,hostname='paddegat')
        
  def get_gravity_at_point(self,radius,x,y,z):
    self.instance.initialize_code()
    self.instance.parameters.epsilon_squared = self.epsilon_squared
    self.instance.particles.add_particles(self.system.particles)
    ax,ay,az=self.instance.get_gravity_at_point(radius,x,y,z)
    self.instance.cleanup_code()
    return ax,ay,az

  def get_potential_at_point(self,radius,x,y,z):
    self.instance.initialize_code()
    self.instance.parameters.epsilon_squared = self.epsilon_squared
    self.instance.particles.add_particles(self.system.particles)
    phi=self.instance.get_potential_at_point(radius,x,y,z)
    self.instance.cleanup_code()
    return phi


def energy_plot(time,ek,ep,eth):
  if not HAS_MATPLOTLIB:
    return
    
  pyplot.figure(figsize = (5, 5))
  pyplot.xlabel(r'time')
  pyplot.ylabel(r'energy')
  pyplot.plot(time,ek)
  pyplot.plot(time,ep)
  pyplot.plot(time,eth)
  pyplot.plot(time,ek+ep+eth)
  test_results_path = get_path_to_results()
  test_results_path = './'
  pyplot.savefig(os.path.join(test_results_path, "evrard_test.png"))

def run_evrard(x):
  conv = nbody_system.nbody_to_si(1000. | units.MSun, 1. | units.parsec)
  evrard=MakeEvrardTest(x,convert_nbody=conv)
  parts=evrard.result
  parts.radius=0. | units.parsec
  parts.h_smooth=0. | units.parsec

  sph=Fi(conv,channel_type='sockets')
  sph.initialize_code()

  sph.parameters.use_hydro_flag=True
  sph.parameters.radiation_flag=False
  sph.parameters.self_gravity_flag=False
  sph.parameters.integrate_entropy_flag=False
  sph.parameters.timestep=conv.to_si(0.0625 | nbody_system.time)  
  sph.parameters.verbosity = 0

  sph.commit_parameters()

  sph.gas_particles.add_particles(parts)
   
  grav=copycat2(Fi, sph, conv, eps2=0.| nbody_system.length**2)
  
  fast=FAST(verbose=False)
  fast.add_system(sph,(grav,),False)
    
  dt=0.0625 | nbody_system.time*0.999999999
  tstart=conv.to_nbody(fast.model_time)
   
  fast.synchronize_model()
  time,Ek,Ep,Eth=[],[],[],[]
  time.append(conv.to_nbody(tstart).value_in(nbody_system.time))
  e=fast.kinetic_energy
  Ek.append(conv.to_nbody(e).value_in(nbody_system.energy))
  e=fast.potential_energy
  Ep.append(conv.to_nbody(e).value_in(nbody_system.energy))
  e=fast.thermal_energy
  Eth.append(conv.to_nbody(e).value_in(nbody_system.energy))
  tnow=tstart
  while tnow<2.99 | nbody_system.time:
    tnow=tnow+dt
    fast.evolve_model(conv.to_si(tnow))
    fast.synchronize_model()
    tcur=fast.model_time
    time.append(conv.to_nbody(tcur).value_in(nbody_system.time))
    print tnow,conv.to_nbody(tcur),conv.to_nbody(sph.model_time)
    e=fast.kinetic_energy
    Ek.append(conv.to_nbody(e).value_in(nbody_system.energy))
    e=fast.potential_energy
    Ep.append(conv.to_nbody(e).value_in(nbody_system.energy))
    e=fast.thermal_energy
    Eth.append(conv.to_nbody(e).value_in(nbody_system.energy))


  time=numpy.array(time)
  Ek=numpy.array(Ek)
  Ep=numpy.array(Ep)
  Eth=numpy.array(Eth)
  energy_plot(time,Ek,Ep,Eth)

if __name__=="__main__":
#  cProfile.run( 'run_evrard(1000)','pprof')  
  run_evrard(10000)
