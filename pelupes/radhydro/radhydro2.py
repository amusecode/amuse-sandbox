import numpy

from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.units import constants

from amuse.community.simplex.interface import SimpleXSplitSet
from amuse.community.sphray.interface import SPHRay

from amuse.community.fi.interface import Fi

from radhydro_code import fi_radhydro

from amuse.datamodel import Particles

from amuse.support.io import write_set_to_file 

suggested_parameter_set=dict()

suggested_parameter_set[Fi]=dict( radiation_flag=False,
                                  self_gravity_flag=False,
                                  gamma=1,
                                  isothermal_flag=True,
                                  integrate_entropy_flag=False )
suggested_parameter_set[SimpleXSplitSet]=dict(number_of_freq_bins=5,
                                              thermal_evolution_flag=1,
                                              blackbody_spectrum_flag=1)
suggested_parameter_set[SPHRay]=dict(default_spectral_type=1,
                                     number_of_rays=1000 | units.Myr**-1)

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
#  gas.flux = 0. | (units.s**-1)
  gas.xion = 0. | units.none

# set source in the center of the box
  sources=Particles(1)
  sources.x=0.*L
  sources.y=0.*L
  sources.z=0.*L
  sources.luminosity=source_luminosity

  return gas,sources

def write_combined_set(filename,radhydro):
    write_set_to_file(radhydro.radhydro_particles_copy(),filename,"amuse", append_to_file=False)    

def T_from_u(xe,u):
    mu=1.| units.amu
    gamma=5./3.    
    return (gamma-1)*mu/(1.+xe)/constants.kB*u

def radhydro_evolve(radhydro,tend,dt,write_snapshots=True):
  t=radhydro.model_time  
  if write_snapshots:
    write_combined_set("radhydro-%6.6i"%0,radhydro) 
  while t<tend-dt/2:
    hydro.evolve_model(t)
    t+=dt

    T=T_from_u(radhydro.rad_particles.xion,radhydro.gas_particles.u)
    print t.in_(units.Myr),
    print "min T:", T.amin().value_in(units.K),
    print "max T:", T.amax().value_in(units.K)

    i=t/dt
    if int(i)%10==0:
      if write_snapshots:
        write_combined_set("radhydro-%6.6i"%int(i),radhydro)

def main( tend=30 | units.Myr,
          N=10000, 
          dt=1. | units.Myr,
          L=15 | units.kpc,    
          source_luminosity=5.e48 | units.s**-1,
          rhoinit=0.001 | (units.amu / units.cm**3),
          Tinit=100. | units.K,
          radhydro_code=fi_radhydro,
          rad_code=SPHRay,
          write_snapshots=True):

  gas,src=iliev_test_5(N=N,L=L,source_luminosity=source_luminosity,
                       rhoinit=rhoinit,Tinit=Tinit)

  converter=nbody_system.nbody_to_si(L**3*rhoinit, L)
  radhydro=radhydro_code(converter, mode='periodic', rad=rad_code)
  hydro_parameters=suggested_parameter_set[hydro_code]
  hydro_parameters['timestep']=dt/2  
  hydro_parameters['periodic_box_size']=2*L
  for x in hydro_parameters:
    radhydro.hydro_parameters.__setattr__(x,hydro_parameters[x])
  radhydro.gas_particles.add_particles(gas)
    
  rad_parameters=suggested_parameter_set[rad_code]
  rad_parameters["box_size"]=2*L
  for x in rad_parameters:
    radhydro.rad_parameters.__setattr__(x,rad_parameters[x])
  radhydro.rad_particles.add_particles(gas)
  radhydro.src_particles.add_particles(src)

  radhydro_evolve(radhydro,tend,dt,write_snapshots=write_snapshots)

if __name__=="__main__":
    from radhydro import main
    main(rad_code=SPHRay) # rad_code=SimpleXSplitSet
