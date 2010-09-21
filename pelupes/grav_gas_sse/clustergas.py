import os
import sys
import numpy

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from grav_gas_sse import grav_gas_sse

from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.units import constants

from amuse.legacy.fi.interface import Fi
from amuse.legacy.gadget2.interface import Gadget2
from amuse.legacy.bhtree.interface import BHTree
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from amuse.legacy.sse.interface import SSE

from amuse.ext.salpeter import SalpeterIMF
from amuse.ext.plummer import MakePlummerModel
from amuse.ext.gasplummer import MakePlummerGasModel
from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.evrard_test import body_centered_grid_unit_cube

import cProfile

numpy.random.seed(123456)

def smaller_nbody_power_of_two(dt, conv):
  nbdt=conv.to_nbody(dt).value_in(nbody_system.time)
  idt=numpy.floor(numpy.log2(nbdt))
  return conv.to_si( 2**idt | nbody_system.time)

def clustergas(sfeff=0.05,Nstar=1000,Ngas=1000, t_end=0.1 | units.Myr,
                 dt_plot= 0.01 | units.Myr, Rscale= 0.3 | units.parsec):

  print "Nstar, Ngas:", Nstar,Ngas
  eps=0.001 * Rscale

  total_star_mass, star_masses = SalpeterIMF().next_set(Nstar)
  total_mass=total_star_mass/sfeff
  print "maxmass:", max(star_masses)
  conv = nbody_system.nbody_to_si(total_mass,Rscale)

  star_parts=MakePlummerModel(Nstar,convert_nbody=conv).result
  star_parts.mass=star_masses
  star_parts.radius=eps

  print "total cluster mass:", total_mass.in_(units.MSun)
  print "star mass:", total_star_mass.in_(units.MSun)
  print "gas mass:", (total_mass-total_star_mass).in_(units.MSun)

  print "t_end:", conv.to_nbody(t_end)
  
  gas_parts=MakePlummerGasModel(Ngas,convert_nbody=conv, base_grid=regular_grid_unit_cube).result
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
  dt_star=dt/4
  dt_sph=dt_star
  dt_fast=dt_star
  dt_feedback=dt
   
  if not dt_star<=dt_fast<=dt_feedback:
    raise Exception 
   
  print 'dt_plot:', conv.to_nbody(dt_plot)
  print 'dt:', conv.to_nbody(dt)
  print 'dt_feedback:', conv.to_nbody(dt_feedback)
  print 'dt_fast:', conv.to_nbody(dt_fast)
  print 'dt_star:', conv.to_nbody(dt_star)
  print 'dt_sph:', conv.to_nbody(dt_sph)
     
  sys=grav_gas_sse(Fi,Fi,SSE,Fi, 
               conv,mgas,star_parts,gas_parts,eps,dt_feedback,dt_fast,
               grav_parameters=(["timestep", dt_star],),
               gas_parameters=(["use_hydro_flag",True],
                               ["radiation_flag",False],
                               ["self_gravity_flag",False],
                               ["verbosity", 0],
                               ["timestep", dt_sph],
                               ["pboxsize", 100 | units.parsec],
#                               ["square_root_timestep_flag",True],
#                               ["sqrt_timestep_crit_constant",.1],
#                               ["acc_timestep_flag",False],
#                               ["gadget_cell_opening_constant",0.01],
#                               ["gadget_cell_opening_flag",True],
                               ["epsilon_squared", eps**2]),
               feedback_efficiency=0.001)

  sys.synchronize_model()
  t=sys.model_time
  tout=t.value_in(units.Myr)
  ek=sys.kinetic_energy.value_in(1.e51*units.erg)
  ep=sys.potential_energy.value_in(1.e51*units.erg)
  eth=sys.thermal_energy.value_in(1.e51*units.erg)
  ef=sys.feedback_energy.value_in(1.e51*units.erg)
  print 't Ek Ep Eth Ef:', tout,ek,ep,eth,ef,ek+ep+eth-ef

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
  
if __name__=="__main__":
  import time
#  cProfile.run("clustergas()","prof")
  t1=time.time()
  clustergas(Ngas=1000)
  t2=time.time()
  print 'time:',  t2-t1
