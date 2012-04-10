import numpy 

from amuse.lab import *
#from amuse.support.units import nbody_system
#from amuse.support.units import units

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

def plot_galaxy(figure, figname, galaxy, bulge, nucleus):
    tl = figure.add_subplot(2,2,1)
    tr = figure.add_subplot(2,2,2)
    bl = figure.add_subplot(2,2,3)
    br = figure.add_subplot(2,2,4)
    tl.set_title('Entire system')
    tr.set_title('Galaxy')
    bl.set_title('stellar cluster')
    br.set_title('black hole sub system')

    com_gal = galaxy.particles.center_of_mass()
    com_bul = bulge.particles.center_of_mass()
    com_nuc = nucleus.particles.center_of_mass()

    tl.plot(galaxy.particles.x.value_in(units.parsec), galaxy.particles.y.value_in(units.parsec), 'g.')
    tl.plot(bulge.particles.x.value_in(units.parsec), bulge.particles.y.value_in(units.parsec), 'b.')
    tl.plot(nucleus.particles.x.value_in(units.parsec), nucleus.particles.y.value_in(units.parsec), 'r.')
    tr.plot(galaxy.particles.x.value_in(units.parsec), galaxy.particles.y.value_in(units.parsec), 'g.')
    bl.plot((bulge.particles.x-com_bul.x).value_in(units.parsec), (bulge.particles.y-com_bul.y).value_in(units.parsec), 'b.')
    bl.plot((nucleus.particles.x-com_bul.x).value_in(units.parsec), (nucleus.particles.y-com_bul.y).value_in(units.parsec), 'r.')
    bl.set_xlim(-1, 1)
    bl.set_ylim(-1, 1)
#    br.plot((nucleus.particles.x-com_nsc.x).value_in(units.parsec), (nucleus.particles.y-com_nsc.y).value_in(units.parsec), 'r.')
    br.plot((nucleus.particles.x).value_in(units.parsec), (nucleus.particles.y).value_in(units.parsec), 'r.')
    br.plot((bulge.particles.x).value_in(units.parsec), (bulge.particles.y).value_in(units.parsec), 'b.')
#    br.set_xlim(-0.002, 0.002)
#    br.set_ylim(-0.002, 0.002)

    figure.savefig(figname)
    figure.clear()

def validate_bridge():

  convert_nsc = nbody_system.nbody_to_si(1000 | units.MSun, 1.0 | units.parsec)
  convert_bul = nbody_system.nbody_to_si(1000 | units.MSun, 1.0 | units.parsec)
  convert_gal = nbody_system.nbody_to_si(1000 | units.MSun, 1.0 | units.parsec)

  total_sim_time= 1 | units.Myr
  print "Total time=", total_sim_time
  sim_clu = convert_bul.to_si(1./16. | nbody_system.time)
  sim_gal = convert_gal.to_si(1./16. | nbody_system.time)
  dt_diag= 2*sim_gal 
  print "time steps:", sim_clu.in_(units.Myr), sim_gal.in_(units.Myr), dt_diag.in_(units.Myr)

  rstep = 1./4. | nbody_system.length
  eps_nuc= convert_nsc.to_si(rstep)
  eps_bul= convert_bul.to_si(rstep)
  eps_gal= convert_gal.to_si(rstep) 

  galaxy=parts_from_file("KingW7N2k.stoa")
  galaxy.radius=0.| nbody_system.length
  galaxy.r = galaxy.position.lengths()
  inner = galaxy.select_array(lambda r:r<0.05|nbody_system.length, ['r',])
  outer = galaxy.select_array(lambda r:r>0.1|nbody_system.length, ['r',])
  middle = galaxy.select(lambda r:r>=0.05|nbody_system.length and r<=0.1|nbody_system.length, ['r',])
  print "galaxy:", len(inner), len(middle), len(outer)
  galaxy = sys_from_parts(BHTree,outer,convert_gal,eps_gal)
#  galaxy = sys_from_parts(PhiGRAPE,outer,convert_gal,eps_gal)
#  bulge = sys_from_parts(PhiGRAPE,middle,convert_bul,eps_bul)
  bulge = sys_from_parts(BHTree,middle,convert_bul,eps_bul)
  nucleus = sys_from_parts(Mikkola,inner,convert_nsc,eps_nuc)
#  nucleus = sys_from_parts(PhiGRAPE,inner,convert_nsc,eps_nuc)

  nuclear_system=bridge.bridge(timestep=sim_clu)
  nuclear_system.add_system(bulge, (nucleus,))
  nuclear_system.add_system(nucleus, (bulge,))

# Second bridge between NSC+cluster and galaxy
  galaxy_system=bridge.bridge(timestep=sim_gal)
  galaxy_system.add_system(galaxy, (nuclear_system,))
  galaxy_system.add_system(nuclear_system, (galaxy,))

  t=0. | units.Myr
  i = 0
  figure = plt.figure(figsize=(6,6))
  E1 = galaxy_system.potential_energy
  while( t < total_sim_time):
    t=t+dt_diag
    galaxy_system.evolve_model(t)
    print_diagnostics(galaxy_system)
    E0 = E1
    E1 = galaxy_system.kinetic_energy+galaxy_system.potential_energy
    dE = E0-E1
    print "Time=", t.in_(units.Myr), "dE=", dE.in_(units.erg), dE/E0
    i+= 1
    figname = "fig_%6.6i.png" %i
    plot_galaxy(figure, figname, galaxy, bulge, nucleus)
  
if __name__ == '__main__':
  validate_bridge()
