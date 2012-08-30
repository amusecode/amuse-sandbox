""" Run test 1 from http://arxiv.org/abs/astro-ph/0603199.  
Note that the box is setup with the origin at the center.  """


import numpy
import argparse


from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.units import constants
from amuse.community.simplex.interface import SimpleXSplitSet
from amuse.community.sphray.interface import SPHRay
from amuse.community.fi.interface import Fi
from amuse.datamodel import Particles
from amuse.support.io import write_set_to_file 

from amuse.ext.evrard_test import regular_grid_unit_cube



import plot_iliev1


OUTPUT_TIME_TOL = 1.0e-4

parser = argparse.ArgumentParser(description='Run test Iliev1.')

parser.add_argument( '-Ngas', 
                     type=int, 
                     default=32768,
                     help='number of gas particles to use' )

parser.add_argument( '-Nngb', 
                     type=int, 
                     default=48,
                     help='number of SPH neighbors to use for hsml' )

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

print
print 'output times: ...'
print OUT_TIMES
print



suggested_parameter_set = {}
Nray_Myr = args.Nray_Myr | units.Myr**-1
suggested_parameter_set[SPHRay] = {'default_spectral_type': -1.,
                                   'number_of_rays': Nray_Myr,
                                   'isothermal_flag': 1 }




def set_gas_and_src( Ngas, Lbox, Lsrc, nHinit, Tinit ):

  """ Given the number of particles and box parameters, sets up gas 
  particles and a source particle. """

  gamma = 5./3.
  mu = 1.0 | units.amu

  gas = Particles(Ngas)

# set particles homogeneously in space from -Lbox/2 to Lbox/2

  if args.grid == 1:

    x,y,z = regular_grid_unit_cube(Ngas).make_xyz()

    gas.x = x * Lbox/2
    gas.y = y * Lbox/2
    gas.z = z * Lbox/2

  elif args.grid == 0:

    gas.x = Lbox/2 * numpy.random.uniform(-1.,1.,Ngas)
    gas.y = Lbox/2 * numpy.random.uniform(-1.,1.,Ngas)
    gas.z = Lbox/2 * numpy.random.uniform(-1.,1.,Ngas)

  else:
    
    print 'args.grid not recognized'
    print 'args.grid = ', args.grid
    sys.exit(1)
  
  

# set other properties

  Nngb = args.Nngb
  Vperpar = Lbox**3 / Ngas * Nngb 
  gas.h_smooth = ( 3 * Vperpar / (4.0 * numpy.pi ) )**(1./3)
  
  gas.u = 1 / (gamma-1) * constants.kB * Tinit/mu

  gas.rho = nHinit

  gas.mass = nHinit*Lbox**3/Ngas

  gas.xion = 1.2e-3

# set source in the center of the box

  sources = Particles(1)
  sources.x = 0.0*Lbox
  sources.y = 0.0*Lbox
  sources.z = 0.0*Lbox
  sources.luminosity = Lsrc

  return gas,sources




def write_rad_set(filename, radpart):
    write_set_to_file(radpart, filename, "amuse", append_to_file=False)    


def rad_evolve( rad, tend, dt, write_snapshots=True ):


  print 'Entered rad evolve'

  # get current model time
  #--------------------------------------------------------
  t = rad.model_time  

  # write initial state if we like
  #--------------------------------------------------------
  if write_snapshots:
    write_rad_set("output/iliev1-%7.7i"%0, rad.gas_particles) 

  # evolve until we get to tend
  #--------------------------------------------------------
  while t < tend-dt/2:
    rad.evolve_model(t+dt)
    t+=dt
    plot_iliev1.quick_check_data( t, rad.gas_particles.mass,
                                  rad.gas_particles.xion )

    plot_iliev1.profile_plots( rad.gas_particles.mass,
                               rad.gas_particles.xion )

    # check if we've reached an output time
    #--------------------------------------------------------
    if write_snapshots:

      t_now = t.value_in(units.Myr)
      t_lbl = t.value_in(units.Myr)*1000
      fname = "output/iliev1-%7.7i"%int(t_lbl)

      for t_check in OUT_TIMES:
        t_diff = numpy.abs( (t_now - t_check) ) 

        if t_diff < OUTPUT_TIME_TOL:

          print fname
          print
          write_rad_set(fname, rad.gas_particles)





def main( Ngas = args.Ngas, 
          tend = args.tend | units.Myr,          
          dt = args.dt | units.Myr,
          Lbox = 13.2 | units.kpc,    
          Lsrc = 5.e48 | units.s**-1,
          nHinit = 1.0e-3 | (units.amu / units.cm**3),
          Tinit = 1.0e4 | units.K,
          rad_code = SPHRay,
          write_snapshots = True ):


  # report to screen
  #---------------------------------------------------------------------
  print 
  print 'Iliev 06 Test 1'
  print 'Ngas: ', Ngas
  print 'tend: ', tend
  print 'dt: ', dt
  print 'Lbox: ', Lbox
  print 'Lsrc: ', Lsrc
  print 'nH: ', nHinit
  print 'Tinit: ', Tinit
  print 'rad code: ', rad_code
  print 


  # set up gas and source particles
  #---------------------------------------------------------------------
  (gas,src) = set_gas_and_src( Ngas, Lbox, Lsrc, nHinit, Tinit )


  # initialize radiative transfer code class
  #---------------------------------------------------------------------
  rad = rad_code( )
  rad_parameters = suggested_parameter_set[rad_code]
  rad_parameters["box_size"]=Lbox

  for x in rad_parameters:
    rad.parameters.__setattr__(x,rad_parameters[x])

  # tell the radiative transfer code about the gas and source particles
  #---------------------------------------------------------------------
  rad.gas_particles.add_particles(gas)
  rad.src_particles.add_particles(src)

  # evolve the system 
  #---------------------------------------------------------------------
  rad_evolve( rad, tend, dt, 
              write_snapshots = write_snapshots )


if __name__=="__main__":
    #from run_iliev1 import main
    main(rad_code=SPHRay) # rad_code=SimpleXSplitSet
