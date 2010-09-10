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

from amuse.legacy.fi.interface import Fi
from amuse.legacy.bhtree.interface import BHTree
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from amuse.legacy.sse.interface import SSE

from amuse.ext.salpeter import SalpeterIMF
from amuse.ext.plummer import MakePlummerModel
from amuse.ext.gasplummer import MakePlummerGasModel
from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.evrard_test import body_centered_grid_unit_cube

import cProfile

def smaller_nbody_power_of_two(dt, conv):
  nbdt=conv.to_nbody(dt).value_in(nbody_system.time)
  idt=numpy.floor(numpy.log2(nbdt))
  return conv.to_si( 2**idt | nbody_system.time)

def clustergas(sfeff=0.1,Nstar=2000,Ngas=100000,  
                 t_end=0.05 | units.Myr,dt_plot= 0.05 | units.Myr):

  total_star_mass, star_masses = SalpeterIMF().next_set(Nstar)
  total_mass=total_star_mass/sfeff

  print "maxmass:", max(star_masses)
  conv = nbody_system.nbody_to_si(total_mass, 1. | units.parsec)

  print "total cluster mass:", total_mass.in_(units.MSun)
  print "star mass:", total_star_mass.in_(units.MSun)
  print "gas mass:", (total_mass-total_star_mass).in_(units.MSun)

  print "t_end:", conv.to_si(t_end).in_(units.Myr)
  
  gas_parts=MakePlummerGasModel(Ngas,convert_nbody=conv, base_grid=regular_grid_unit_cube).result
  gas_parts.h_smooth=0. | units.parsec
  gas_parts.mass=gas_parts.mass*(1-sfeff)

  print "gas particle mass:",  ((total_mass-total_star_mass)/len(gas_parts)).in_(units.MSun)

  mgas=(total_mass-total_star_mass)/len(gas_parts)
  print max(gas_parts.u)**0.5
  star_parts=MakePlummerModel(Nstar,convert_nbody=conv).result
  star_parts.radius=0. | units.parsec
  star_parts.mass=star_masses

  eps=0.001| units.parsec
  star_parts.radius=eps
  
  dt =smaller_nbody_power_of_two(dt_plot, conv)
  dt_star=dt/8
  dt_sph=dt_star
  dt_fast=4*dt_star
  dt_feedback=dt_fast
   
  print 'dt_plot:', conv.to_nbody(dt_plot)
  print 'dt:', conv.to_nbody(dt)
  print 'dt_feedback:', conv.to_nbody(dt_feedback)
  print 'dt_fast:', conv.to_nbody(dt_fast)
  print 'dt_star:', conv.to_nbody(dt_star)
  print 'dt_sph:', conv.to_nbody(dt_sph)
     
  sys=grav_gas_sse(PhiGRAPE,Fi,SSE,Fi, 
               conv,mgas,star_parts,gas_parts,eps,dt_feedback,dt_fast,
               grav_parameters=(),
               gas_parameters=(["use_hydro_flag",True],
                               ["radiation_flag",False],
                               ["self_gravity_flag",False],
                               ["verbosity", 0],
                               ["timestep", dt_sph]),
               feedback_efficiency=0.01)

  t=sys.model_time
  while (t<t_end):
    sys.evolve_model(t+dt)
    t=sys.model_time
  
if __name__=="__main__":
  cProfile.run("clustergas()","prof")
#  clustergas()
