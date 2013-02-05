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

from radiativehydro import RadiativeHydro

from amuse.datamodel import Particles
from amuse.support.io import write_set_to_file 

from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.evrard_test import body_centered_grid_unit_cube


import plot_iliev5



OUTPUT_TIME_TOL = 1.0e-4
numpy.random.seed(1234567)

CC = plot_iliev5.Iliev5Vars()

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

parser.add_argument( '-grid', 
                     type=int, 
                     default=1,
                     help='particles set to ( 1=grid, 0=random )' )


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
#                               'isothermal_flag':False,
                               'integrate_entropy_flag':True,
                               'periodic_box_size': CC.Lbox}



# SPHray
#-------------------------------------------------------------------
Nray_Myr = args.Nray_Myr | units.Myr**-1
suggested_parameter_set[SPHRay] = {'default_spectral_type': 1,
                                   'number_of_rays': Nray_Myr,
                                   'isothermal_flag': 0,
                                   'box_size': CC.Lbox,
                                   'boundary_condition': 0}


# SimpleXSplitSet
#-------------------------------------------------------------------
suggested_parameter_set[SimpleXSplitSet] = {'number_of_freq_bins':5,
                                            'thermal_evolution_flag':1,
                                            'blackbody_spectrum_flag':1,
                                            'metal_cooling_flag':0,
                                            'box_size': CC.Lbox }



def set_gas_and_src( Ngas, Lbox, Lsrc, rho_init, T_init ):

  """ Given the number of particles and box parameters, sets up gas 
  particles and a source particle. """

  gamma = 5./3.
  mu = 1.0 | units.amu

# set particles homogeneously in space from -Lbox/2 to Lbox/2
# NOTE reset Ngas

  if args.grid == 1:

    x,y,z = regular_grid_unit_cube(Ngas).make_xyz()
    #x,y,z = body_centered_grid_unit_cube(Ngas).make_xyz()

    Ngas = len(x)
    gas = Particles(Ngas)

    gas.x = x * Lbox/2
    gas.y = y * Lbox/2
    gas.z = z * Lbox/2


  elif args.grid == 0:

    gas = Particles(Ngas)

    gas.x = Lbox/2 * numpy.random.uniform(-1.,1.,Ngas)
    gas.y = Lbox/2 * numpy.random.uniform(-1.,1.,Ngas)
    gas.z = Lbox/2 * numpy.random.uniform(-1.,1.,Ngas)

  else:
    
    print 'args.grid not recognized'
    print 'args.grid = ', args.grid
    sys.exit(1)



# set zero velocities

  gas.vx = 0. | (units.km/units.s)
  gas.vy = 0. | (units.km/units.s)
  gas.vz = 0. | (units.km/units.s)

# set other properties

  gas.h_smooth = 0.0 * Lbox # will be set in call to hydro code
  
  gas.u = 1 / (gamma-1) * constants.kB * T_init / mu

  gas.rho = rho_init

  gas.mass = rho_init*Lbox**3/Ngas

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
          Lbox = 30.0 | units.kpc,    
          Lsrc = 5.e48 | units.s**-1,
          rho_init = 1.0e-3 | (units.amu / units.cm**3),
          T_init = 1.0e2 | units.K,
          rad_code = SPHRay,
          #rad_code = SimpleXSplitSet,
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
  print 'rho_init: ', rho_init
  print 'T_init: ', T_init
  print 'rad code: ', rad_code
  print 'hydro code: ', hydro_code
  print 


  # set up gas and source particles
  #---------------------------------------------------------------------
  (gas,src) = set_gas_and_src( Ngas, Lbox, Lsrc, rho_init, T_init )

  print 'gas x min/max: ', gas.x.min().value_in( units.kpc ), \
      gas.x.max().value_in( units.kpc )

  print 'gas y min/max: ', gas.y.min().value_in( units.kpc ), \
      gas.y.max().value_in( units.kpc )

  print 'gas z min/max: ', gas.z.min().value_in( units.kpc ), \
      gas.z.max().value_in( units.kpc )

  print 'gas vx min/max: ', gas.vx.min().value_in( units.km/units.s ), \
      gas.vx.max().value_in( units.km/units.s )

  print 'gas vy min/max: ', gas.vy.min().value_in( units.km/units.s ), \
      gas.vy.max().value_in( units.km/units.s )

  print 'gas vz min/max: ', gas.vz.min().value_in( units.km/units.s ), \
      gas.vz.max().value_in( units.km/units.s )



  print


  # set up a system of units in which G=1.  the function takes any two 
  # fundamental dimensions (in this case mass and length) and calculates
  # the third such that G = 1. 
  #---------------------------------------------------------------------
  converter = nbody_system.nbody_to_si( Lbox**3 * rho_init, Lbox )


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

  for x in hydro_parameters:
    radhydro.hydro_parameters.__setattr__(x,hydro_parameters[x])

  radhydro.gas_particles.add_particles(gas)

  # set rad code parameters
  #---------------------------------------------------------------------
  rad_parameters = suggested_parameter_set[rad_code]
  if rad_code == SPHRay:
    rad_parameters["momentum_kicks_flag"] = False

  for x in rad_parameters:
    radhydro.rad_parameters.__setattr__(x,rad_parameters[x])

  radhydro.rad_particles.add_particles(gas)
  radhydro.src_particles.add_particles(src)


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

    print t.in_(units.Myr)
    print "min/max T:", T.amin().in_(units.K), T.amax().in_(units.K)
    print "min/max x_ion:", xion.min(), xion.max()
    print "min/max v:", v.amin().in_(units.kms), v.amax().in_(units.kms)
    print

    # check if we've reached an output time
    #--------------------------------------------------------
    if write_snapshots:

      t_now = t.value_in(units.Myr)
      t_lbl = t.value_in(units.Myr)*1000
      fname = "output/iliev5-%7.7i"%int(t_lbl)

      for t_check in OUT_TIMES:
        t_diff = numpy.abs( (t_now - t_check) ) 

        if t_diff < OUTPUT_TIME_TOL:

          # function combines data from hydro particles
          # and radiation particles
          
          gas_particles = radhydro.radhydro_particles_copy()

          write_set_to_file( gas_particles,
                             fname, "amuse", append_to_file=False )

          plot_iliev5.plot_profiles( data=gas_particles, t=t )
          plot_iliev5.plot_images( data=gas_particles, t=t )





if __name__=="__main__":
    #from run_iliev5 import main
    #main(hydro_code=Fi, rad_code=SPHRay) # rad_code=SimpleXSplitSet
    main(hydro_code=Fi, rad_code=SimpleXSplitSet) # rad_code=SimpleXSplitSet


    # note do not pass uniform grid distribution to Simplex !!!
    # DO pass the boxsize to Simplex
