"""
STONES IN CLUSTER
stars (gravity and stellar evolution) and zero-mass particles (stones)
  -- stars: integrated with N-body (time-step --heta),
            stellar evolution (times-step --tse)
  -- stones: advance in the potential bridged to the stars (timestep --brdt)
"""

import numpy
import sys
import os

from amuse.units import units, constants, nbody_system, quantities
from amuse.units.optparse import OptionParser
from amuse.datamodel import ParticlesSuperset, Particles
from amuse.ic.plummer import new_plummer_sphere
from amuse.community.sse.interface import SSE
#from amuse.community.seba.interface import SeBa
#from amuse.community.hermite0.interface import Hermite
from amuse.community.huayno.interface import Huayno
#from amuse.ic.brokenimf import new_broken_power_law_mass_distribution
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.io import write_set_to_file
from amuse.couple import bridge

# high order bridge
# from amuse.ext.composition_methods import *
"""
  included:
    LEAPFROG
    SPLIT_4TH_S_M6
    SPLIT_4TH_S_M5
    SPLIT_4TH_S_M4
    SPLIT_6TH_SS_M11
    SPLIT_6TH_SS_M13
    SPLIT_8TH_SS_M21
    SPLIT_10TH_SS_M35
"""

def get_cluster(n_stars, 
                radius=2.5|units.parsec,
                vir=0.5,
                m_min=0.3|units.MSun, m_max=25.0|units.MSun):
  """ 
  ini star cluster
    Plummer sphere + IMF
  """
  
  #masses = new_broken_power_law_mass_distribution(n_stars, 
                                                  #mass_boundaries=[0.08, 0.5, 100] | units.MSun, 
                                                  #alphas=[-1.3,-2.3])
  masses = new_salpeter_mass_distribution(n_stars, 
                                          mass_min=m_min,
                                          mass_max=m_max)
  converter = nbody_system.nbody_to_si(masses.sum(),radius)
  stars = new_plummer_sphere(n_stars, converter)
  stars.mass = masses
  stars.move_to_center()
  stars.scale_to_standard(converter, virial_ratio=vir)
  
  #sigma = (converter.to_si(1./numpy.sqrt(2.) | (nbody_system.length/nbody_system.time))).in_(units.kms)
  sigma = numpy.mean(stars.velocity.lengths().value_in(units.kms))
  
  print ' ** cluster: m_tot =', stars.mass.sum().in_(units.MSun), 'm_min =', masses.min().in_(units.MSun), 'm_max =', masses.max().in_(units.MSun), ', <v> =', sigma, 'kms'
  
  return stars, converter

def get_stones(n_stones,
               stars,
               converter,
               radius=2.5|units.parsec,
               vir=0.5):
  """
  ini stones
    Plummer sphere (same as for stars), zero-mass particles
  """
  
  stones = new_plummer_sphere(n_stones, converter)
  stones.mass = 0.0 | units.MJupiter
  stars_aux = stars.copy()
  st_st = ParticlesSuperset([stars_aux, stones])
  st_st.scale_to_standard(converter, virial_ratio=vir)
  
  sigma = numpy.mean(stones.velocity.lengths().value_in(units.kms))
  
  print " ** stones: m_tot =", stones.mass.sum().in_(units.MSun), ', <v> =', sigma, 'kms'
  
  return stones

def grav_code(conv, heta):
  """ ini Huayno for stars """
  code = Huayno(conv, channel_type="sockets", mode="openmp")
  code.parameters.epsilon_squared = (0.0 | units.parsec)**2
  code.parameters.timestep_parameter = heta
  return code

def se_code():
  """ ini SE code for stars """
  code = SSE()
  #code = SeBa()
  return code

def softening(particles):
  """
  --- from Carmen ---
  optimum softening lenght. Formula valid
  only for clusters in virial equlibrium
  Ref: Aarseth book
  """
  N = len(particles.mass)
  U = particles.potential_energy()
  Rvir = 0.5*constants.G*particles.mass.sum()**2/abs(U)
  epsilon = 4.0*Rvir/N
  return epsilon

class advance_without_selfgravity(object):
    """
    to advance particles
    """
    def __init__(self, particles, time= 0 |units.Myr):
        self.particles = particles
        self.model_time = time
    
    def evolve_model(self, t_end):
        dt = t_end - self.model_time
        self.particles.position += self.particles.velocity*dt
        self.model_time= t_end
    
    @property
    def potential_energy(self):
        return quantities.zero
    
    @property 
    def kinetic_energy(self):
        return (0.5*self.particles.mass*self.particles.velocity.lengths()**2).sum()

def evolve_with_stellar_evolution(tend,
                                  converter, 
                                  grav, se,
                                  timestep=0.25|units.Myr,
                                  E0=None):
  """
  evolve system
    -- stellar evolution -- mass is updated every timestep
  """
  
  ### channels
  channel_from_se_to_grav = se.particles.new_channel_to(grav.particles)
  
  ### initial quantities
  if E0 is None:
    E0 = grav.particles.kinetic_energy() + grav.particles.potential_energy()
  
  ### evolve
  t = grav.model_time
  while (t < (tend-0.5*timestep)):
    
    grav.evolve_model(t+timestep)
    t = grav.model_time
    se.evolve_model(t)
    channel_from_se_to_grav.copy_attributes(["mass"])
    
    mass_t= grav.particles.mass.sum()
    E_t = grav.kinetic_energy + grav.potential_energy
    dE = (E_t - E0)/ E0
    
    print "\t", t.in_(units.Myr), mass_t.in_(units.MSun), dE
    
  return t, E0

def print_option_parser(parser_options, file_par):
  """ print options into file """
  f_par = open(file_par, 'w')
  for oi in parser_options.__dict__:
    #print oi, parser_options.__dict__[oi]
    f_par.write(oi+' '+str(parser_options.__dict__[oi])+'\n')
  f_par.close()
  return
  
    
def new_option_parser():
  result = OptionParser()
  result.add_option("--n1", 
                    dest="n_stars", type="int", default = 10,
                    help="number of stars in the cluster [%default]")
  result.add_option("--n2", 
                    dest="n_stones", type="int", default = 10,
                    help="number of stones in the cluster [%default]")
  result.add_option("--tend", 
                    dest="tend", default = 1|units.Myr, type="float", unit=units.Myr,
                    help="evolution time [%default]")
  result.add_option("--tsnap", 
                    dest="tsnap", default = 1|units.Myr, type="float", unit=units.Myr,
                    help="snapshot time-step [%default]")
  result.add_option("--tse", 
                    dest="tse", default = 0.05|units.Myr, type="float", unit=units.Myr,
                    help="stellar evolution time-step [%default]")
  result.add_option("--brdt", 
                    dest="brdt", default = 0.001|units.Myr, type="float", unit=units.Myr,
                    help="bridge time-step to advance stones [%default]")
  result.add_option("--fout", 
                    dest="fout", default="snap",
                    help="snapshot file name [%default]")
  result.add_option("--dsnap", 
                    dest="dsnap", default="./snap/",
                    help="snapshot directory [%default]")
  result.add_option("--rc", 
                    dest="r_stars", default=1.0|units.parsec, type="float", unit=units.parsec,
                    help="cluster radius [%default]")
  result.add_option("--heta", 
                    dest="heta", default=0.01, type="float",
                    help="Huayno time-step parameter [%default]")
  result.add_option("--vir", 
                    dest="vir", default=0.5, type="float",
                    help="virial ratio [%default]")
  result.add_option("--seed", 
                    dest="seed", default=666, type="int",
                    help="set seed for the IC [%default]")
  result.add_option("--m_min", 
                    dest="m_min", default = 0.3|units.MSun, type="float", unit=units.MSun,
                    help="minimal mass [%default]")
  result.add_option("--m_max", 
                    dest="m_max", default = 25.0|units.MSun, type="float", unit=units.MSun,
                    help="maxima mass [%default]")
  return result

if __name__ in ('__main__', '__plot__'):
  o, arguments  = new_option_parser().parse_args()
  
  print_option_parser(o, 'stones_in_cluster.par')
  
  if o.seed is not None:
    numpy.random.seed(seed=o.seed)
  
  
  if (o.tse < o.brdt):
    print " ** warning : SE times-step (tse) must be at least the time-step to advance stones (brdt) ! ** "
    sys.exit()
  else:
    pass
  
  try:
     os.mkdir(o.dsnap)
  except:
     print " ** ", o.dsnap, "already exists"
  
  ### IC for the star cluster
  stars, converter = get_cluster(o.n_stars, radius=o.r_stars, 
                                 vir=o.vir, 
                                 m_min=o.m_min, m_max=o.m_max)
  stars.collection_attributes.timestamp = 0.|units.Myr
  
  ### initialize gravity for star cluster
  grav = grav_code(converter, o.heta)
  #soft = softening(stars)
  #grav.parameters.epsilon_squared = soft**2
  grav.particles.add_particles(stars)
  grav.commit_particles()
  #print grav.get_timestep_parameter()
  
  ### initialize stellar evolution for the cluster
  se = se_code()
  se.particles.add_particles(stars)
  se.commit_particles()
  
  ### IC for stones
  stones = get_stones(o.n_stones, 
                      stars,
                      converter, radius=o.r_stars,
                      vir=o.vir)
  stones.collection_attributes.timestamp = 0.|units.Myr
  
  ### initialize stones under the influence of the cluster
  advance_stones = advance_without_selfgravity(stones)
  cluster_and_stones = bridge.Bridge(timestep=o.brdt)
  cluster_and_stones.add_system(advance_stones, (grav,), False)    # cluster acts on the stones
  cluster_and_stones.add_system(grav, (), False)                   # cluster evolves on its own
  ###
  
  #
  # just in case of high-order bridge is needed
  # cluster_and_stones = bridge.Bridge(timestep=o.brdt, method=SPLIT_4TH_S_M6)
  #
  
  ### snapshots timesteps
  snap_times = (numpy.append( numpy.arange(0., o.tend.value_in(units.Myr), o.tsnap.value_in(units.Myr)), \
                o.tend.value_in(units.Myr)) ) | units.Myr
  print " ** snapshot times:", snap_times
  
  ### channels
  # to update star cluster
  channel_from_grav_to_stars = grav.particles.new_channel_to(stars)
  # to update stones
  channel_from_advance_to_stones = advance_stones.particles.new_channel_to(stones)
  
  ### stars and stones together
  bodies = ParticlesSuperset([stars, stones])
  
  i = 0
  E0 = None
  for snap_time in snap_times:
    
     t_ev, E0 = evolve_with_stellar_evolution(snap_time,
                                          converter, 
                                          cluster_and_stones, 
                                          se,
                                          timestep=o.tse,
                                          E0=E0)
     
     channel_from_grav_to_stars.copy()
     channel_from_advance_to_stones.copy()
     
     bodies.collection_attributes.timestamp = t_ev
     
     fout_grav = o.dsnap+o.fout+"_gr_{0:04d}.hdf5".format(i)
     fout_se   = o.dsnap+o.fout+"_se_{0:04d}.hdf5".format(i)
     print "\t", snap_time, " -> ", fout_grav
     write_set_to_file(bodies, fout_grav,"amuse")
     write_set_to_file(se.particles, fout_se, "amuse")
     i += 1
  
  grav.stop()
  se.stop()

  ### ___ it's over now ___
