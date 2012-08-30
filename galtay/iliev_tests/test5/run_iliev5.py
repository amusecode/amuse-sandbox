""" Run test 5 from http://arxiv.org/abs/0905.2920.  
Note that the box is setup with the origin at the center.  """


import numpy
import argparse


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

OUTPUT_TIME_TOL = 1.0e-4
numpy.random.seed(1234567)


parser = argparse.ArgumentParser(description='Run test Iliev5.')

parser.add_argument( '-Ngas', 
                     type=int, 
                     default=32768,
                     help='number of gas particles to use' )

parser.add_argument( '-Nray_Myr', 
                     type=int, 
                     default=1000,
                     help='number of rays to trace / Myr' )

parser.add_argument( '-tend', 
                     type=float, 
                     default=500.0,
                     help='time to run test for [Myr]' )

parser.add_argument( '-dt', 
                     type=float, 
                     default=1.0,
                     help='time between sucessive rad evolve calls [Myr]' )

parser.add_argument( '-tout', 
                     type=float, 
                     default=[10.0, 30.0, 100.0, 200.0, 500.0],
                     nargs='+',
                     help='list of output times [Myr]' )

parser.add_argument( '-movie_dt', 
                     type=float, 
                     default=None,
                     help='make output every MOVIE_DT Myr for making movies' )


args = parser.parse_args()

if args.movie_dt:
  nouts = numpy.int( args.tend / args.movie_dt ) + 1
  OUT_TIMES = numpy.linspace( 0.0, args.tend, nouts )
  OUT_TIMES = OUT_TIMES[1:]
else:
  OUT_TIMES = args.tout


print OUT_TIMES



# set up default parameters for different codes
#======================================================================
suggested_parameter_set = {}


# Gadget2
#-------------------------------------------------------------------
suggested_parameter_set[Gadget2] = {}

# Fi
#-------------------------------------------------------------------
suggested_parameter_set[Fi] = {'radiation_flag':False,
                               'self_gravity_flag':False,
#                               'gamma':1,
#                               'isothermal_flag':True,
                               'integrate_entropy_flag':True}

# SPHray
#-------------------------------------------------------------------
Nray_Myr = args.Nray_Myr | units.Myr**-1
suggested_parameter_set[SPHRay] = {'default_spectral_type': 1,
                                   'number_of_rays': Nray_Myr,
                                   'isothermal_flag': 0 }


# SimpleXSplitSet
#-------------------------------------------------------------------
suggested_parameter_set[SimpleXSplitSet] = {'number_of_freq_bins':5,
                                            'thermal_evolution_flag':1,
                                            'blackbody_spectrum_flag':1,
                                            'metal_cooling_flag':0 }



def set_gas_and_src( Ngas, Lbox, Lsrc, nHinit, Tinit ):

  """ Given the number of particles and box parameters, sets up gas 
  particles and a source particle. """

  gamma = 5./3.
  mu = 1.0 | units.amu

  gas = Particles(Ngas)

# set particles homogeneously in space from -Lbox/2 to Lbox/2

  gas.x = Lbox/2 * numpy.random.uniform(-1.,1.,Ngas)
  gas.y = Lbox/2 * numpy.random.uniform(-1.,1.,Ngas)
  gas.z = Lbox/2 * numpy.random.uniform(-1.,1.,Ngas)

# set zero velocities

  gas.vx = 0. | (units.km/units.s)
  gas.vy = 0. | (units.km/units.s)
  gas.vz = 0. | (units.km/units.s)

# set other properties

  gas.h_smooth = 0.0 * Lbox # will be set in call to hydro code
  
  gas.u = 1 / (gamma-1) * constants.kB * Tinit / mu

  gas.rho = nHinit

  gas.mass = nHinit*Lbox**3/Ngas

  gas.xion = 1.2e-3

  #gas.dudt=gas.u/(1.| units.Myr)

# set source in the center of the box

  sources = Particles(1)
  sources.x = 0.0*Lbox
  sources.y = 0.0*Lbox
  sources.z = 0.0*Lbox
  sources.luminosity = Lsrc

  return gas,sources




def T_from_u(xe, u):
  mu=1.| units.amu
  gamma = 5./3.    
  T = (gamma-1)*mu / (1.+xe) / constants.kB * u
  return T




def main( Ngas = args.Ngas, 
          tend = args.tend | units.Myr,          
          dt = args.dt | units.Myr,
          Lbox = 15.0 | units.kpc,    
          Lsrc = 5.e48 | units.s**-1,
          nHinit = 1.0e-3 | (units.amu / units.cm**3),
          Tinit = 1.0e2 | units.K,
          rad_code = SPHRay,
          hydro_code = Fi,
          write_snapshots = True ):


  # report to screen
  #---------------------------------------------------------------------
  print 
  print 'Iliev 09 Test 5'
  print 'Ngas: ', Ngas
  print 'tend: ', tend
  print 'dt: ', dt
  print 'Lbox: ', Lbox
  print 'Lsrc: ', Lsrc
  print 'nH: ', nHinit
  print 'Tinit: ', Tinit
  print 'rad code: ', rad_code
  print 'hydro code: ', hydro_code
  print 


  # set up gas and source particles
  #---------------------------------------------------------------------
  (gas,src) = set_gas_and_src( Ngas, Lbox, Lsrc, nHinit, Tinit )


  # set up a system of units in which G=1.  the function takes any two 
  # fundamental dimensions (in this case mass and length) and calculates
  # the third such that G = 1. 
  #---------------------------------------------------------------------
  converter = nbody_system.nbody_to_si( Lbox**3 * nHinit, Lbox )


  # initialize the radhydro class
  #---------------------------------------------------------------------
  def hydrocode():
    return hydro_code(converter, mode='periodic')
  
  def radcode():
    return rad_code()

  radhydro = RadiativeHydro( rad=radcode, hydro=hydrocode )

  # set hydro code parameters
  #---------------------------------------------------------------------
  hydro_parameters = suggested_parameter_set[hydro_code]
  hydro_parameters['timestep'] = dt/2  
  hydro_parameters['periodic_box_size'] = Lbox

  for x in hydro_parameters:
    radhydro.hydro_parameters.__setattr__(x,hydro_parameters[x])

  radhydro.gas_particles.add_particles(gas)

  # set rad code parameters
  #---------------------------------------------------------------------
  rad_parameters = suggested_parameter_set[rad_code]
  rad_parameters["box_size"] = Lbox
  rad_parameters["momentum_kicks_flag"] = False

  for x in rad_parameters:
    radhydro.rad_parameters.__setattr__(x,rad_parameters[x])

  radhydro.rad_particles.add_particles(gas)


  # evolve system
  #---------------------------------------------------------------------
  t = radhydro.model_time  
  while t < tend-dt/2:
    t += dt
    radhydro.evolve_model(t)

    T = T_from_u(radhydro.rad_particles.xion, radhydro.gas_particles.u)
    xion = radhydro.rad_particles.xion
    v = (radhydro.gas_particles.vx**2 + 
         radhydro.gas_particles.vz**2 + 
         radhydro.gas_particles.vy**2 )**0.5

    print t.in_(units.Myr),
    print "min T:", T.amin().in_(units.K),
    print "max T:", T.amax().in_(units.K),
    print "min x_ion:", xion.min(),
    print "max x_ion:", xion.max(),
    print "min v:", v.amin().in_(units.kms),
    print "max v:", v.amax().in_(units.kms)
    

  write_set_to_file( radhydro.radhydro_particles_copy(),
                     'gas_final',"amuse", append_to_file=False )


if __name__=="__main__":
    #from run_iliev1 import main
    main(hydro_code=Fi, rad_code=SPHRay) # rad_code=SimpleXSplitSet
