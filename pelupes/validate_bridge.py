import numpy 

from amuse.support.units import nbody_system
from amuse.support.units import units

from amuse.support.data import core 
from amuse.support.io import text

from amuse.legacy.hermite0.interface import Hermite
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from amuse.legacy.bhtree.interface import BHTree
from amuse.legacy.fi.interface import Fi
import bridge

from amuse.support.codes import channel

channel.MessageChannel.DEBUGGER=channel.MessageChannel.XTERM

def parts_from_file(input_file):
  io=text.TableFormattedText(input_file)
  io.header_prefix_string=";"
  io.attribute_names=["ignore","mass","x","y","z","vx","vy","vz"]
  io.attribute_types=[units.none,nbody_system.mass,nbody_system.length,
    nbody_system.length,nbody_system.length,nbody_system.speed,
    nbody_system.speed,nbody_system.speed]
  return io.load()

def sys_from_parts(base_class,parts,converter,eps=None,timestep=None,usegl=False):
  if usegl:
    interface=base_class(converter,use_gl=True)
  else:
    interface=base_class(converter)
  interface.initialize_code()
#  if interface.__class__.__name__=="Fi":
#    value=converter.to_si(nbody_system.length).in_(units.kpc).number 
#    interface.set_unitl_in_kpc(value)
#    value=converter.to_si(nbody_system.mass).in_(units.MSun).number 
#    interface.set_unitm_in_msun(value)    
  if eps is not None:
    interface.parameters.epsilon_squared = eps**2 
  if timestep is not None:
    interface.parameters.timestep = timestep 
  interface.particles.add_particles(parts)
  return interface

def sys_from_file(base_class,input_file,converter, eps=None,timestep=None,usegl=False):
  parts=parts_from_file(input_file)
  parts.radius=0.| nbody_system.length
  return sys_from_parts(base_class,parts,converter,eps,timestep,usegl)

def shift_sys(system,dx,dy,dz,dvx,dvy,dvz):
  parts=system.particles.copy()
  parts.x=parts.x+dx
  parts.y=parts.y+dy
  parts.z=parts.z+dz
  parts.vx=parts.vx+dvx
  parts.vy=parts.vy+dvy
  parts.vz=parts.vz+dvz
  parts.copy_values_of_state_attributes_to(system.particles)

def validate_bridge():
 
  convert_clu = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
  convert_gal = nbody_system.nbody_to_si(1.e7 | units.MSun, 10. | units.parsec)

  convert_gal_to_clu=lambda x: convert_clu.to_nbody(convert_gal.to_si(x))
  convert_clu_to_gal=lambda x: convert_gal.to_nbody(convert_clu.to_si(x))

  total_sim_time=convert_clu_to_gal(200 | nbody_system.time)
  print "Total time=", total_sim_time
  sim_step=1./128 | nbody_system.time
  dt_diag=1 | nbody_system.time  

  eps_clu= convert_gal_to_clu(2.e-4 | nbody_system.length)
  eps_gal= 6.25e-3 | nbody_system.length

  cluster=sys_from_file(PhiGRAPE, "king7_2k.txt", convert_clu, eps_clu,usegl=True)
  galaxy=sys_from_file(Fi,"king9_100k.txt",convert_gal,eps_gal,sim_step,usegl=True)

  shift_sys(cluster,
    convert_clu.to_si( convert_gal_to_clu(2.5 | nbody_system.length) ),0 | units.m,0 | units.m,
    0| units.m/units.s,convert_clu.to_si( convert_gal_to_clu(0.65 | nbody_system.speed) ),0| units.m/units.s )
  
  if hasattr(cluster,"start_viewer"): cluster.start_viewer()
  if hasattr(galaxy,"viewer"): galaxy.viewer()

    
  bridgesys=bridge.bridge()
  bridgesys.add_system(galaxy)
  bridgesys.add_system(cluster, (galaxy,) )

  print "hellow",galaxy.parameters.timestep.in_(units.Myr)
  print convert_gal.to_si(nbody_system.time)
  print convert_gal.to_si(nbody_system.time).in_(units.Myr)

  print galaxy.get_unitl_in_kpc()
  print galaxy.get_unitm_in_msun()

  t=0. | nbody_system.time
  while( t < total_sim_time):
    t=t+dt_diag
    bridgesys.evolve( convert_gal.to_si(t),
                    timestep=convert_gal.to_si( sim_step )) 
  
if __name__ == '__main__':
  validate_bridge()
