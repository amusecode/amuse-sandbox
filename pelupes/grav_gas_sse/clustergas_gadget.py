import os
import sys
import numpy

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from grav_gas_sse import grav_gas_sse

from amuse.units import nbody_system
from amuse.units import units
from amuse.units import constants

from amuse.community.fi.interface import Fi
from boxedfi import BoxedFi
from amuse.community.octgrav.interface import Octgrav
from amuse.community.gadget2.interface import Gadget2
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPE

from SSEplus import SSEplus

from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.evrard_test import body_centered_grid_unit_cube

import cProfile

#numpy.random.seed(12345)

from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import salpeter_masses
def smaller_nbody_power_of_two(dt, conv):
  nbdt=conv.to_nbody(dt).value_in(nbody_system.time)
  idt=numpy.floor(numpy.log2(nbdt))
  return conv.to_si( 2**idt | nbody_system.time)


def clustergas_restart(runid, snapshot, t_end=30. | units.Myr,
                         dt_plot= 0.05 | units.Myr,newid=None,
                         new_gas_options=()):
  if runid is None:
    raise Exception                       
  print "restarting run "+runid+" at snapshot %i"%snapshot                       
  print
  
  conv,sys=grav_gas_sse.load_system_state(runid+"/dump-%6.6i" %snapshot,
             new_gas_options)

  dt =smaller_nbody_power_of_two(dt_plot, conv)

  if newid is not None:
    runid=newid
    print "switch to:",runid
    try:
      os.mkdir(runid)
    except:
      raise Exception

  print
  print "so far so good...resuming evolve!"

  nsnap=snapshot
  t=sys.model_time
  time_offset=sys.time_offset
  
  while (t < t_end-dt/2):
    sys.evolve_model(t+dt)
    sys.synchronize_model()
    t=sys.model_time
    tout=(t).value_in(units.Myr)
    ek=sys.kinetic_energy.value_in(1.e51*units.erg)
    ep=sys.potential_energy.value_in(1.e51*units.erg)
    eth=sys.thermal_energy.value_in(1.e51*units.erg)
    ef=sys.feedback_energy.value_in(1.e51*units.erg)
    print 't Ek Ep Eth Ef:', tout,ek,ep,eth,ef,ek+ep+eth-ef
    nsnap+=1
    sys.dump_system_state(runid+"/dump-%6.6i" %nsnap)
                         

def clustergas(sfeff=0.05,Nstar=1000,Ngas=1000, t_end=30. | units.Myr,
                 dt_plot= 0.05 | units.Myr, Rscale= 0.5 | units.parsec,
                 runid="runtest",feedback_efficiency=0.01, subvirialfac=1.):

  try:
    os.mkdir(runid)
  except:
    pass  
  print "Nstar, Ngas:", Nstar,Ngas
  eps=0.05 * Rscale
  eps_star=0.001 * Rscale

  star_masses = new_salpeter_mass_distribution(Nstar, mass_max = 100. | units.MSun)
  total_star_mass = star_masses.sum()
  total_mass=total_star_mass/sfeff
  print "maxmass:", max(star_masses)
  print "Tcross", (2.8*(Rscale**3/ constants.G/ total_star_mass)**0.5).in_(units.Myr)

  print sorted(star_masses)[-10:]
#  raise Exception

  conv = nbody_system.nbody_to_si(total_mass,Rscale)

  star_parts = new_plummer_model(Nstar,convert_nbody=conv)
  star_parts.mass=star_masses
  star_parts.radius=eps_star
  
# sub-virialized
  star_parts.vx=star_parts.vx*subvirialfac
  star_parts.vy=star_parts.vy*subvirialfac
  star_parts.vz=star_parts.vz*subvirialfac

  print "total cluster mass:", total_mass.in_(units.MSun)
  print "star mass:", total_star_mass.in_(units.MSun)
  print "gas mass:", (total_mass-total_star_mass).in_(units.MSun)

  print "t_end:", conv.to_nbody(t_end)
  
  gas_parts=new_plummer_gas_model(Ngas,convert_nbody=conv, base_grid=body_centered_grid_unit_cube)
  gas_parts.h_smooth=0. | units.parsec
  gas_parts.mass=gas_parts.mass*(1-sfeff)

  print "gas particle mass:",  ((total_mass-total_star_mass)/len(gas_parts)).in_(units.MSun)

  mu=1.4 | units.amu
  gamma1=1.6667-1
#  print 'min Temp:', (gamma1*min(gas_parts.u)*(1.4*units.amu)/constants.kB).in_(units.K)
  print 'min Temp:', (gamma1*min(gas_parts.u)*mu/constants.kB).in_(units.K)

  mgas=(total_mass-total_star_mass)/len(gas_parts)
  print max(gas_parts.u)**0.5
  
  dt =smaller_nbody_power_of_two(dt_plot, conv)
  dt_star=smaller_nbody_power_of_two( 0.001 | units.Myr, conv) #0.001
  dt_sph=dt_star
  dt_fast=dt_star*8  # 8
  dt_feedback=dt_fast*2 # 1
    
  if not dt_star<=dt_fast<=dt_feedback:
    raise Exception 
     
  print 'dt_plot:', conv.to_nbody(dt_plot), dt_plot.in_(units.Myr)
  print 'dt:', conv.to_nbody(dt), dt.in_(units.Myr)
  print 'dt_feedback:', conv.to_nbody(dt_feedback), dt_feedback.in_(units.Myr)
  print 'dt_fast:', conv.to_nbody(dt_fast), dt_fast.in_(units.Myr)
  print 'dt_star:', conv.to_nbody(dt_star), dt_star.in_(units.Myr)
  print 'dt_sph:', conv.to_nbody(dt_sph), dt_sph.in_(units.Myr)

  star_parts.move_to_center()
  gas_parts.move_to_center()
  
  sys=grav_gas_sse(PhiGRAPE,Gadget2,SSEplus,Octgrav, 
               conv,mgas,star_parts,gas_parts,dt_feedback,dt_fast,
               grav_parameters=(["epsilon_squared", eps_star**2],["timestep_parameter", 0.001]),
#                                ["timestep", dt_star]),
               gas_parameters=(["time_max", 32. | units.Myr],
#                               ["courant", 0.15 | units.none],
                               ["n_smooth", 64 | units.none],
#                               ["artificial_viscosity_alpha",1.| units.none],
                               ["n_smooth_tol", 0.005 |units.none],
                               ## NB
                               ["max_size_timestep",7500. | units.yr],
                               ## NB
                               ["time_limit_cpu", 3600000 | units.s]),
               couple_parameters=(["epsilon_squared", eps**2],
                            ["opening_angle",0.5]),
               feedback_efficiency=feedback_efficiency, feedback_radius=0.025*Rscale)

  nsnap=0
  sys.synchronize_model()
  t=sys.model_time
  tout=t.value_in(units.Myr)
  ek=sys.kinetic_energy.value_in(1.e51*units.erg)
  ep=sys.potential_energy.value_in(1.e51*units.erg)
  eth=sys.thermal_energy.value_in(1.e51*units.erg)
  ef=sys.feedback_energy.value_in(1.e51*units.erg)
  print 't Ek Ep Eth Ef:', tout,ek,ep,eth,ef,ek+ep+eth-ef
  sys.dump_system_state(runid+"/dump-%6.6i" %nsnap)

  print 'max smooth:', max(sys.sph.gas_particles.radius).in_(units.parsec)
  print 'mean smooth:', sys.sph.gas_particles.radius.mean().in_(units.parsec)
  print 'min smooth:', min(sys.sph.gas_particles.radius).in_(units.parsec)
  print 'eps:',eps.in_(units.parsec)
  print 'feedback radius:',(0.01*Rscale).in_(units.parsec)

  while (t<t_end-dt/2):
    sys.evolve_model(t+dt)
    sys.synchronize_model()
    t=sys.model_time
    tout=t.value_in(units.Myr)
    ek=sys.kinetic_energy.value_in(1.e51*units.erg)
    ep=sys.potential_energy.value_in(1.e51*units.erg)
    eth=sys.thermal_energy.value_in(1.e51*units.erg)
    ef=sys.feedback_energy.value_in(1.e51*units.erg)
    print 't Ek Ep Eth Ef:', tout,ek,ep,eth,ef,ek+ep+eth-ef
    nsnap+=1
    sys.dump_system_state(runid+"/dump-%6.6i" %nsnap)
  
if __name__=="__main__":
  import time
#  cProfile.run("clustergas()","prof")
  t1=time.time()
  nb=clustergas(Nstar=1000,Ngas=10000)
  t2=time.time()
  print 'time:',  t2-t1
