import numpy 

from amuse.support.units import nbody_system
from amuse.support.units import units

from amuse.support.data import core 
from amuse.support.io import text

from amuse.community.mikkola.interface import Mikkola
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.bhtree.interface import BHTree
import bridge

from matplotlib import pyplot as plt
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
  if converter == None:
    interface=base_class()
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

def shift_sys(system,dpos,dvel):
  parts=system.particles.copy()
  parts.x=parts.x+dpos[0]
  parts.y=parts.y+dpos[1]
  parts.z=parts.z+dpos[2]
  parts.vx=parts.vx+dvel[0]
  parts.vy=parts.vy+dvel[1]
  parts.vz=parts.vz+dvel[2]
  parts.copy_values_of_all_attributes_to(system.particles)


def print_diagnostics(galaxy_system):
  Ek = galaxy_system.kinetic_energy
  Ep = galaxy_system.potential_energy
  print "Energies (at t=",galaxy_system.time,": ", Ek, Ep, Ek+Ep, "Q=", Ek/Ep

def validate_bridge():

  convert_nsc = nbody_system.nbody_to_si(100 | units.MSun, 0.001 | units.parsec)
  convert_clu = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
  convert_gal = nbody_system.nbody_to_si(1.e7 | units.MSun, 10. | units.parsec)

  total_sim_time= 1 | units.Myr
  print "Total time=", total_sim_time
  sim_clu = convert_clu.to_si(1./16 | nbody_system.time)
  sim_gal = convert_gal.to_si(1./16 | nbody_system.time)
  dt_diag= 2*sim_gal 
  print "time steps:", sim_clu.in_(units.Myr), sim_gal.in_(units.Myr), dt_diag.in_(units.Myr)

  eps_nsc= convert_nsc.to_si(1.0/4. | nbody_system.length)
  eps_clu= convert_clu.to_si(1.0/4. | nbody_system.length)
  eps_gal= convert_gal.to_si(1.0/4. | nbody_system.length)

  nuclear_cl=sys_from_file(PhiGRAPE, "PlummerN10.stoa", convert_nsc, eps_nsc)
  cluster=sys_from_file(PhiGRAPE, "KingW7N2k.stoa", convert_clu, eps_clu)
  galaxy=sys_from_file(BHTree, "KingW7N2k.stoa", convert_gal, eps_gal)

  dpos = nbody_system.length([2.5, 0, 0]) 
  dvel = nbody_system.speed([0, 0.65, 0]) 
  shift_sys(nuclear_cl, convert_clu.to_si(dpos), convert_clu.to_si(dvel))

  nuclear_system=bridge.bridge(timestep=sim_clu)
  nuclear_system.add_system(cluster, (nuclear_cl,))
  nuclear_system.add_system(nuclear_cl, (cluster,))

  dpos = nbody_system.length([1.5, 0, 1.0]) 
  dvel = nbody_system.speed([0, 1.3, 0]) 
  # shift the nuclear system with respect to the galaxy in dpos and dvel
  shift_sys(nuclear_system, convert_gal.to_si(dpos), convert_gal.to_si(dvel))

# Second bridge between NSC+cluster and galaxy
  galaxy_system=bridge.bridge(timestep=sim_gal)
  galaxy_system.add_system(galaxy, (nuclear_system,))
  galaxy_system.add_system(nuclear_system, (galaxy,))

#  galaxy_system.evolve_model()
#  plt.plot(galaxy_system.particles.x.value_in(units.parsec), galaxy_system.particles.y.value_in(units.parsec), 'r.')

  t=0. | units.Myr
  i = 0
  figure = plt.figure(figsize=(6,6))
  while( t < total_sim_time):
    print "Time=", t.in_(units.Myr)
    t=t+dt_diag
    galaxy_system.evolve_model(t)
    print_diagnostics(galaxy_system)

    plot = figure.add_subplot(1,1,1)
    plot.plot(galaxy.particles.x.value_in(units.parsec), galaxy.particles.y.value_in(units.parsec), 'g.')
    plot.plot(cluster.particles.x.value_in(units.parsec), cluster.particles.y.value_in(units.parsec), 'b.')
    plot.plot(nuclear_cl.particles.x.value_in(units.parsec), nuclear_cl.particles.y.value_in(units.parsec), 'r.')
    i+= 1
    figname = "fig_%6.6i.png" %i
    figure.savefig(figname)
    figure.clear()

#  plt.show()


#  print "hellow",galaxy.parameters.timestep.in_(units.Myr)
#  print convert_gal.to_si(nbody_system.time)
#  print convert_gal.to_si(nbody_system.time).in_(units.Myr)

#  print galaxy.get_unitl_in_kpc()
#  print galaxy.get_unitm_in_msun()

  
if __name__ == '__main__':
  validate_bridge()
