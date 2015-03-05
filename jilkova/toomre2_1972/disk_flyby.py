import numpy 
import os
import time

from amuse.units.optparse import OptionParser
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSuperset
from amuse.community.huayno.interface import Huayno
from amuse.community.kepler_orbiters.interface import Kepler
from amuse.io import write_set_to_file, read_set_from_file
from amuse.couple import bridge

from gravity_subsets import GravitySubsets

#from plot_snaps_all_02 import plot_all_snaps

# high order bridge
#from amuse.ext.composition_methods import *

# optimization
#import cProfile

pi_180 = numpy.pi/180.0

def rm_file(file_name):
  """
  delete the file 'file_name', if it exists
  """  
  if(os.path.isfile(file_name)):
    os.system('rm %s' % (file_name))

def initialize_gravity_subsets(bodies, converter, gravity_code, particles_per_subset,
                               **options):
  """
  initialize gravity code to integrate particles in subsets
  """
  subsets_gravity = GravitySubsets(gravity_code, converter, bodies, particles_per_subset, **options)
  subsets_gravity.add_particles()     # perticles used for the initialization added (bodies)
  subsets_gravity.commit_particles()
  
  return subsets_gravity

def integrate_2disks_flyby(stars, planetesimals, t_end, n_steps, 
                           file_out, file_redir, 
                           huayno_eta_nbody,
                           huayno_eps2,
                           huayno_inttype_parameter,
                           nsub,
                           verbose=True):
  """
  Integrate a binary orbit with test particles (e.g., disks around stars of the binary).
  -- integrated with gravity in subsets
     (the planetesimals are followed by N-body code in subsets including the binary)
  """
  
  # to get the total computational time of the script
  t0 = time.time()
  
  # calculate number of steps for each integration
  t_ini = 0.0 | units.yr       # initial time
  t_tot = t_end - t_ini        # total time 
  dt = t_tot / float(n_steps)  # time step
  
  if verbose==True:
    print " ** n_steps: tot = ", n_steps
    
  if file_out is not None:
    print " ** output file exists -- removing", file_out
    rm_file(file_out)
    
  # set the converter
  converter=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
  
  ### n-body integration
  # make superset
  bodies = ParticlesSuperset([stars, planetesimals])
  # initialize HUAYNO for superset
  gravity_nbody = initialize_gravity_subsets(bodies, converter, Huayno, nsub, 
                                             timestep_parameter=huayno_eta_nbody, 
                                             epsilon_squared=huayno_eps2,
                                             inttype_parameter=huayno_inttype_parameter)
  
  if verbose==True:
    t11= time.time()
    dt = t11-t0
    dt_kepler1 = dt
  
  stars, planetesimals, t_stop_nbody = evolve_disk_flyby_together_in_subsets(stars, planetesimals, gravity_nbody, 
                                                                             t_ini, t_end, n_steps, 
                                                                             converter, file_out, verbose=verbose)
  bodies_fin = ParticlesSuperset([stars, planetesimals])
  
  if verbose is True:
    t1= time.time()
    dt = t1-t0
    print "Performace data: N =", len(stars)+len(planetesimals), "dt=", dt, "s =", dt/60.0, "min"
  
  return bodies_fin
  
def evolve_disk_flyby_together_in_subsets(stars, planetesimals, gravity, 
                               t_start, t_end, n_steps, converter, file_out,
                               verbose=True):

  bodies = ParticlesSuperset([stars, planetesimals])
  channel_from_gr_to_framework = gravity.particles.new_channel_to(bodies)
  Etot_init = gravity.kinetic_energy() + gravity.potential_energy()
  Etot = Etot_init
  
  duration = t_end - t_start
  dt = duration / float(n_steps-1)
  time = 0.0 | units.yr

  if verbose is True:
    #print " ** Huayno timestep parameter =", gravity.list_of_instances[0].get_timestep_parameter()
    print " ** evolving: t_start = ", t_start.in_(units.Myr), "t_end = ", t_end.in_(units.Myr), ", dt = ", dt.in_(units.Myr)
    print " \t", "time", "\t\t", "dE"

  while time<=duration:
    gravity.evolve_model(time)
    #channel_from_gr_to_framework.copy()
    
    # this is because the channel does not work
    bodies = gravity.particles.copy()
    
    bodies.collection_attributes.timestamp = time + t_start
    
    Ekin = gravity.kinetic_energy()
    Epot = gravity.potential_energy()
    Etot = Ekin + Epot
    dE = Etot_init-Etot
    nb_E = converter.to_nbody(Etot)
    if verbose is True:
      print " \t\t", time+t_start, "\t", dE/Etot_init#, "\t", nb_E
    
    if file_out is not None:
      write_set_to_file(bodies, file_out, "hdf5")
    
    time += dt
  
  gravity.stop()
  
  # this is because the channel does not work
  stars = bodies[:2].copy()
  planetesimals = bodies[2:].copy()
  
  return stars, planetesimals, bodies.collection_attributes.timestamp
