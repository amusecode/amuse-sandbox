import os
import sys
import numpy

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    
from amuse.support.units import nbody_system
from amuse.support.units import units

    
from amuse.community.fi import interface as interface
from amuse.community.gadget2.interface import Gadget2
from amuse.ext.evrard_test import MakeEvrardModel
from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.evrard_test import body_centered_grid_unit_cube
from amuse.test.amusetest import get_path_to_results

import logging
#logging.basicConfig(level=logging.DEBUG)

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
  pyplot.xlim(0,3)
  test_results_path = get_path_to_results()
  test_results_path = './'
  pyplot.savefig(os.path.join(test_results_path, "evrard_test.png"))

def run_evrard(x):
  convert_nbody = nbody_system.nbody_to_si(1000. | units.MSun, 1. | units.parsec)
  evrard=MakeEvrardModel(x,convert_nbody=convert_nbody)
  parts=evrard.result
  parts.h_smooth=0. | units.parsec

  sph=interface.Fi(convert_nbody,use_gl=False)
#  sph=Gadget2(convert_nbody)
  sph.initialize_code()

  sph.parameters.use_hydro_flag=True
  sph.parameters.radiation_flag=False
  sph.parameters.integrate_entropy_flag=True
  sph.parameters.artificial_viscosity_alpha=0.5
  sph.parameters.courant=0.15
  sph.parameters.timestep=0.05 | nbody_system.time  

  sph.commit_parameters()

  sph.gas_particles.add_particles(parts)
    
  if hasattr(sph,"viewer"):
    sph.viewer()

  dt=0.0625 | nbody_system.time
  tstart=convert_nbody.to_nbody(sph.model_time) 
  sph.synchronize_model()

  time,Ek,Ep,Eth=[],[],[],[]
  time.append(convert_nbody.to_nbody(tstart).value_in(nbody_system.time))
  e=sph.kinetic_energy
  Ek.append(convert_nbody.to_nbody(e).value_in(nbody_system.energy))
  e=sph.potential_energy
  Ep.append(convert_nbody.to_nbody(e).value_in(nbody_system.energy))
  e=sph.thermal_energy
  Eth.append(convert_nbody.to_nbody(e).value_in(nbody_system.energy))
  tnow=tstart
  while tnow<2.99 | nbody_system.time:
    tnow=tnow+dt
    sph.evolve_model(tnow)
    sph.synchronize_model()
    tcur=sph.model_time
    time.append(convert_nbody.to_nbody(tcur).value_in(nbody_system.time))
    e=sph.kinetic_energy
    Ek.append(convert_nbody.to_nbody(e).value_in(nbody_system.energy))
    e=sph.potential_energy
    Ep.append(convert_nbody.to_nbody(e).value_in(nbody_system.energy))
    e=sph.thermal_energy
    Eth.append(convert_nbody.to_nbody(e).value_in(nbody_system.energy))


  time=numpy.array(time)
  Ek=numpy.array(Ek)
  Ep=numpy.array(Ep)
  Eth=numpy.array(Eth)
  energy_plot(time,Ek,Ep,Eth)

if __name__=="__main__":
  run_evrard(1000)
