import numpy

from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.units import constants

from amuse.community.simplex.interface import SimpleXSplitSet
from amuse.community.sphray.interface import SPHRay

from amuse.community.fi.interface import Fi

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
  gas.flux = 0. | (units.s**-1)
  gas.xion = 0. | units.none

# set source in the center of the box
  sources=Particles(1)
  sources.x=0.*L
  sources.y=0.*L
  sources.z=0.*L
  sources.luminosity=source_luminosity

  return gas,sources

def write_combined_set(filename,hydropart,radpart):
    p=radpart.copy()
    channel = hydropart.new_channel_to(p)
    channel.copy_attributes(["mass","x","y","z","vx","vy","vz","rho","u"])  	  
    write_set_to_file(p,filename,"amuse", append_to_file=False)    

def T_from_u(xe,u):
    mu=1.| units.amu
    gamma=5./3.    
    return (gamma-1)*mu/(1.+xe)/constants.kB*u

def radhydro_evolve(hydro,rad,tend,dt,write_snapshots=True):
  t=rad.model_time  
  update_rad_from_hydro_channel = hydro.gas_particles.new_channel_to(rad.gas_particles)
  update_hydro_from_rad_channel = rad.gas_particles.new_channel_to(hydro.gas_particles)
  if write_snapshots:
    write_combined_set("radhydro-%6.6i"%0,hydro.gas_particles,rad.gas_particles) 
  while t<tend-dt/2:
    hydro.evolve_model(t+dt/2)    
    update_rad_from_hydro_channel.copy_attributes(["x","y","z","rho","h_smooth"])
    rad.evolve_model(t+dt)
    update_hydro_from_rad_channel.copy_attributes(["u"])
    hydro.evolve_model(t+dt)    
    t+=dt

    T=T_from_u(rad.gas_particles.xion,rad.gas_particles.u)
    print t.in_(units.Myr),        
    print "min T:", T.amin().value_in(units.K),
    print "max T:", T.amax().value_in(units.K)

    i=t/dt
    if int(i)%10==0:
      if write_snapshots:
        write_combined_set("radhydro-%6.6i"%int(i),hydro.gas_particles,rad.gas_particles)

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
  hydro=hydro_code(converter, mode='periodic')
  hydro_parameters=suggested_parameter_set[hydro_code]
  hydro_parameters['timestep']=dt/2  
  hydro_parameters['periodic_box_size']=2*L
  for x in hydro_parameters:
    hydro.parameters.__setattr__(x,hydro_parameters[x])
  hydro.gas_particles.add_particles(gas)
    
  rad=rad_code()
  rad_parameters=suggested_parameter_set[rad_code]
  rad_parameters["box_size"]=2*L
  for x in rad_parameters:
    rad.parameters.__setattr__(x,rad_parameters[x])
  rad.gas_particles.add_particles(gas)
  rad.src_particles.add_particles(src)

  radhydro_evolve(hydro,rad,tend,dt,write_snapshots=write_snapshots)

if __name__=="__main__":
    from radhydro import main
    main (rad_code=SimpleXSplitSet)  # (rad_code=SPHRay)
