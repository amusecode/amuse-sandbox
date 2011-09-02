import numpy 

from amuse.support.data import core 
from amuse.units import nbody_system
from amuse.units import units
from amuse.io import text
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.bhtree.interface import BHTree
import bridge

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

def test():
  convert = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
  eps=2.e-4 | nbody_system.length

  test_class=BHTree

  cluster=sys_from_file(test_class, "king7_2k.txt", convert, eps)
  cluster.synchronize_model()
  print convert.to_nbody(cluster.potential_energy)
  print convert.to_nbody(cluster.kinetic_energy)

  parts=cluster.particles.copy()
  parts1=parts.select_array(lambda x: (x > 0 | units.m), ['x'] )
  parts2=parts.select_array(lambda x: (x < 0 | units.m), ['x'] )
  cluster1=sys_from_parts(test_class, parts1, convert, eps)
  cluster2=sys_from_parts(test_class, parts2, convert, eps)

  print
  cluster1.synchronize_model()
  print convert.to_nbody(cluster1.potential_energy)
  print convert.to_nbody(cluster1.kinetic_energy)
  cluster2.synchronize_model()
  print convert.to_nbody(cluster2.potential_energy)
  print convert.to_nbody(cluster2.kinetic_energy)
  print
  
  bridgesys=bridge.bridge()
  bridgesys.add_system(cluster1, (cluster2,) )
  bridgesys.add_system(cluster2, (cluster1,) )

  print convert.to_nbody(bridgesys.potential_energy)
  print convert.to_nbody(bridgesys.kinetic_energy)
  
if __name__=="__main__":
  test()  
