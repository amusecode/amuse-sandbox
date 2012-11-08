import os.path
import numpy
from time import sleep

import cPickle as pickle

from matplotlib import pyplot

from amuse.plot import loglog, xlabel, ylabel
from amuse.io import write_set_to_file, read_set_from_file
from amuse.community.mesa.interface import MESA
from amuse.community.fi.interface import Fi
from amuse.units import units, generic_unit_system, nbody_system, constants
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.datamodel import Particle, Particles, ParticlesSuperset


from amuse.ext.star_to_sph import convert_stellar_model_to_SPH

from amuse.ext.sph_to_star import SPH2StellarModel, convert_SPH_to_stellar_model

def new_particles(N=1000,M=1.0| units.MSun,age=5.|units.Gyr):
  
  input_file = "star_"+str(N)+'_'+str(M.value_in(units.MSun))+'_'+str(age.value_in(units.Gyr))
  if os.path.exists(input_file):
      return read_set_from_file(input_file,"amuse")
  
  stellar_evolution = MESA()
  stellar_evolution.particles.add_particle(Particle(mass=M))
  stellar_evolution.evolve_model(age)
  particles = convert_stellar_model_to_SPH(
      stellar_evolution.particles[0],
      N, 
      seed=12345, base_grid_options=dict(type="sobol")
  ).gas_particles
  stellar_evolution.stop()
  
  hydrodynamics = Fi(nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.RSun))
  hydrodynamics.parameters.timestep = 100.0 | units.s
#  hydrodynamics = Gadget2(ConvertBetweenGenericAndSiUnits(1.0|units.MSun, 1.0|units.RSun, 1.0e3|units.s))
  hydrodynamics.gas_particles.add_particles(particles)
  print "evolving"
  hydrodynamics.evolve_model(1.0|units.s)
  hydrodynamics.gas_particles.copy_values_of_attributes_to(["density", "u", "pressure"], particles)
  hydrodynamics.stop()
  write_set_to_file(particles, input_file,"amuse")
  return particles

def hr_tracks(stellar_evolution):
    stars=stellar_evolution.particles
    data=dict()
    
    i=0
    for star in stars:
        stardata=data.setdefault(i,dict())
        i+=1
        stardata['luminosity'] = [] | units.LSun
        stardata['temperature'] = [] | units.K
        stardata['time'] = [] | units.yr
        star.evolve_for(star.age+ (50. | units.Myr))
        while star.age < 14. | units.Gyr:
#        while star.stellar_type < 10 | units.stellar_type:
#        while star.age < 10 | units.Gyr:
#        while star.luminosity < 10.| units.LSun:
            stardata['luminosity'].append(star.luminosity)
            stardata['temperature'].append(star.temperature)
            stardata['time'].append(star.age)
            star.evolve_one_step()
#            star.evolve_for(star.time_step*5)
            print star.stellar_type, star.age, star.mass, star.luminosity, star.radius        
    return data
    
def plot_track(data):
    colors=["r","g","b","c","k"]
    pyplot.figure(figsize = (8, 6))
#    pyplot.title('Hertzsprung-Russell diagram', fontsize=12)
    for star,stardata in data.items():
        temperature=stardata['temperature']
        luminosity=stardata['luminosity']
        loglog(temperature, luminosity,colors[star])
        if stardata['time'][0] > 100 | units.Myr:
          loglog(temperature[0:1],luminosity[0:1],'x'+colors[star])
        else:
          loglog(temperature[0:1],luminosity[0:1],'^'+colors[star])  
#        text(1.25*temperature[-1], 0.5*luminosity[-1], mass)
    xlabel('Effective Temperature')
    ylabel('Luminosity')
    pyplot.xlim(10**6, 10**3)
    pyplot.ylim(1.0e-1,1.e4)
    pyplot.savefig("mesasph_hr2_100k.eps")   
    
if __name__ in ('__main__', '__plot__'):
    
  age=5. | units.Gyr  
  N=100000
    
  stellar_evolution = MESA()
#  stellar_evolution.parameters.stabilize_new_stellar_model_flag=False
  stellar_evolution.particles.add_particle(Particle(mass=1.0|units.MSun)) # reference particle

  for age in [5.|units.Gyr,8.|units.Gyr,10.|units.Gyr]:#,10000,10000]:
    parts=new_particles(N,age=age)
    model = convert_SPH_to_stellar_model(parts) # model is from center to surface
    stellar_evolution.new_particle_from_model(model, age)
  data=hr_tracks(stellar_evolution)
  plot_track(data)

  f=open('data_100k_5_8_10','wb')
  pickle.dump(data,f)
  f.close()
