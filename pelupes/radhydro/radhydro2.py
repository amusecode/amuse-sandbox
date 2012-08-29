import numpy

numpy.random.seed(1234567)

from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.units import constants

from amuse.community.simplex.interface import SimpleXSplitSet
from amuse.community.sphray.interface import SPHRay

from amuse.community.fi.interface import Fi
from amuse.community.gadget2.interface import Gadget2

from radhydro_code import RadiativeHydro

from amuse.datamodel import Particles

from amuse.support.io import write_set_to_file 

suggested_parameter_set=dict()

suggested_parameter_set[Fi]=dict( radiation_flag=False,
                                  self_gravity_flag=False,
                                  gamma=1,
                                  isothermal_flag=True,
                                  integrate_entropy_flag=False )

suggested_parameter_set[Gadget2]=dict( )

suggested_parameter_set[SimpleXSplitSet]=dict(number_of_freq_bins=5,
                                              thermal_evolution_flag=1,
                                              blackbody_spectrum_flag=1)
suggested_parameter_set[SPHRay]=dict(default_spectral_type=1,
                                     number_of_rays=1000 | units.Myr**-1,
                                     isothermal_flag=False)

def iliev_test_5( N=10000,
                  L=15. | units.kpc,
                  source_luminosity=5.e48 | units.s**-1,
                  rhoinit=0.001 | (units.amu / units.cm**3),
                  Tinit=100. | units.K,
                  hydro_parameters=dict(),
                  rad_parameters=dict()):

  gamma=5./3.
  mu=1.| units.amu

  gas=Particles(N)
# set particles homogeneously in space
  gas.x=L*numpy.random.uniform(-1.,1.,N)
  gas.y=L*numpy.random.uniform(-1.,1.,N)
  gas.z=L*numpy.random.uniform(-1.,1.,N)
# set other properties
  gas.h_smooth = 0. | units.parsec
  gas.vx = 0. | (units.km/units.s)
  gas.vy = 0. | (units.km/units.s)
  gas.vz = 0. | (units.km/units.s)
  gas.u = 1/(gamma-1)*constants.kB * Tinit/mu
  gas.rho = rhoinit
  gas.mass = rhoinit*(2*L)**3/N
  gas.xion = 0.1 | units.none
  gas.dudt=gas.u/(1.| units.Myr)

# set source in the center of the box
  sources=Particles(1)
  sources.x=0.*L
  sources.y=0.*L
  sources.z=0.*L
  sources.luminosity=source_luminosity

  return gas,sources

def T_from_u(xe,u):
    mu=1.| units.amu
    gamma=5./3.    
    return (gamma-1)*mu/(1.+xe)/constants.kB*u

def main( tend=30 | units.Myr,
          N=10000, 
          dt=1. | units.Myr,
          L=15 | units.kpc,    
          source_luminosity=5.e48 | units.s**-1,
          rhoinit=0.001 | (units.amu / units.cm**3),
          Tinit=100. | units.K,
          hydro_code=Fi,
          rad_code=SPHRay,
          write_snapshots=True):

  gas,src=iliev_test_5(N=N,L=L,source_luminosity=source_luminosity,
                       rhoinit=rhoinit,Tinit=Tinit)

  converter=nbody_system.nbody_to_si(L**3*rhoinit, L)
  def hydrocode():
    return hydro_code(converter, mode='periodic')
  
  def radcode():
    return rad_code()
  
  radhydro=RadiativeHydro(rad=radcode,hydro=hydrocode)
  
  hydro_parameters=suggested_parameter_set[hydro_code]
  hydro_parameters['timestep']=dt/2  
  hydro_parameters['periodic_box_size']=2*L
  for x in hydro_parameters:
    radhydro.hydro_parameters.__setattr__(x,hydro_parameters[x])
  radhydro.gas_particles.add_particles(gas)
      
  rad_parameters=suggested_parameter_set[rad_code]
  rad_parameters["box_size"]=2*L
  rad_parameters["momentum_kicks_flag"]=False
  for x in rad_parameters:
    radhydro.rad_parameters.__setattr__(x,rad_parameters[x])
  radhydro.rad_particles.add_particles(gas)
#  radhydro.src_particles.add_particles(src)

  t=radhydro.model_time  
  while t<tend-dt/2:
    t+=dt
    radhydro.evolve_model(t)

    T=T_from_u(radhydro.rad_particles.xion,radhydro.gas_particles.u)
    xion=radhydro.rad_particles.xion
    v=(radhydro.gas_particles.vx**2+radhydro.gas_particles.vz**2+radhydro.gas_particles.vy**2)**0.5
    print t.in_(units.Myr),
    print "min T:", T.amin().in_(units.K),
    print "max T:", T.amax().in_(units.K),
    print "min x_ion:", xion.min(),
    print "max x_ion:", xion.max(),
    print "min v:", v.amin().in_(units.kms),
    print "max v:", v.amax().in_(units.kms)
    
  write_set_to_file(radhydro.radhydro_particles_copy(),'gas_final',"amuse", append_to_file=False)
    
if __name__=="__main__":
    main(hydro_code=Fi,rad_code=SPHRay,Tinit=10000. | units.K) # rad_code=SimpleXSplitSet
  
