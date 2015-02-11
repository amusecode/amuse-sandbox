import numpy

from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSuperset
from amuse.support.options import option, OptionalAttributes
from amuse.units import units

class GravitySubsets():
  
  def __init__(self,
               gravity_code,
               converter,
               binary_and_stones,
               particles_per_subset,
               binary_diff_tol=1.e-7,   # should correspond to the desired energy conservation
               **options):
    
    self.gravity_code = gravity_code
    self.converter = converter
    self.binary_and_stones = binary_and_stones
    self.particles_per_subset = particles_per_subset
    self.binary_diff_tol = binary_diff_tol
    self.options = options
    
    self.number_of_subsets = int(numpy.ceil((len(binary_and_stones) - 2.0)/particles_per_subset))
    self.model_time = 0.0 | units.Myr
    self.binary, self.list_binary_and_stones_subsets = self.divide_to_subsets()
    
    self.particles = binary_and_stones.copy()
    
    self.list_of_instances = self.initialize_code()
    #print self.list_of_instances
    #for instance_i in self.list_of_instances:
      #print instance_i.get_inttype_parameter()
    
    # testing the subset division
    #print self.binary
    #for subset_i in self.list_binary_and_stones_subsets:
      #print len(subset_i), subset_i.total_mass().in_(units.MSun)
  
  def divide_to_subsets(self):
    """
    divide the set to subsets
      each subset includes the first two particles, wchich are assumed to be the 
        two most massive ones
      ( each subsets includes the two most massive particles and )
      particles_per_subset of the remaining particles
      (the last subset can have less particles)
    """
    #print len(self.binary_and_stones)
    #binary_and_stones_sorted_by_mass = self.binary_and_stones.sorted_by_attribute('mass')
    #binary = binary_and_stones_sorted_by_mass[-2:]
    binary = self.binary_and_stones[:2]
    stones = self.binary_and_stones - binary
    #print len(self.binary_and_stones)
    binary_and_stones_subsets = []
    stones_last = stones
    for i in range(self.number_of_subsets-1):
      i0 = i * self.particles_per_subset
      i1 = (i+1) * self.particles_per_subset
      stones_and_binary_i = binary + stones[i0:i1]
      binary_and_stones_subsets.append(stones_and_binary_i)
      stones_last = stones_last - stones[i0:i1]
      #print len(stones_last)
    binary_and_stones_subsets.append(binary + stones_last)
    return binary, binary_and_stones_subsets
  
  def run_instances_in_serial(self, end_time):
    # final particle set of stones
    stones_final = Particles(0)
    
    # list of final binary particles sets
    #  should be the same binary number_of_subsets-times
    list_binary_final = []
    
    for i,instance_i in enumerate(self.list_of_instances):
      instance_i.evolve_model(end_time)
      stones_final.add_particles(instance_i.particles[2:].copy())
      list_binary_final.append(instance_i.particles[:2].copy())
    
    # check if binary evolved the same in all subsets
    diff_flag = self.final_binary_diff(list_binary_final)
    if diff_flag!=0:
      print "\t binary of the first subset used for the evolution"
    
    # final particle set after evolve
    binary_and_stones_final = Particles(0)
    binary_and_stones_final.add_particles(list_binary_final[0])
    binary_and_stones_final.add_particles(stones_final)
    
    return binary_and_stones_final
  
  def final_binary_diff(self, list_binary_final):
    """
    check if binary position and velocity evolved the same in all subsets
      calculate relative (with respect to the lengths of the position, resp. velocity, 
        vector of the first binary in the list) difference in each component of 
        position and velocity vectors for each subset and compare these with 
        binary_diff_tol parameter
      warning given in case the difference is higher
      returns diff_flag: 0 -- all differences < binary_diff_tol
        1 -- difference in position > binary_diff_tol
        2 -- difference in velocity > binary_diff_tol
        3 -- difference in position and velocity > binary_diff_tol
    """
    diff_flag = 0
    
    delta_pos = []
    delta_vel = []
    
    for bin_i in list_binary_final[1:]:
      #print (bin_i.position-list_binary_final[0].position).in_(units.AU)
      #print (bin_i.velocity-list_binary_final[0].velocity).in_(units.kms)
      delta_pos.append((bin_i[0].position-list_binary_final[0][0].position)/list_binary_final[0][0].position.lengths())
      delta_vel.append((bin_i[0].velocity-list_binary_final[0][0].velocity)/list_binary_final[0][0].velocity.lengths())
      delta_pos.append((bin_i[1].position-list_binary_final[0][1].position)/list_binary_final[0][1].position.lengths())
      delta_vel.append((bin_i[1].velocity-list_binary_final[0][1].velocity)/list_binary_final[0][1].velocity.lengths())
    
    delta_pos_m = numpy.array(delta_pos)
    delta_pos_m = delta_pos_m[~numpy.isnan(delta_pos_m)]
    delta_pos_over_tol = (abs(delta_pos_m)>self.binary_diff_tol)
    if (delta_pos_over_tol.any() == True):
      print " !!! warning: relative position difference after gravity_subset higher then", self.binary_diff_tol, "!!! "
      print "\t binary position vector differs by:"
      print "\t", delta_pos_m[delta_pos_over_tol]
      diff_flag = 1
    
    delta_vel_m = numpy.array(delta_vel)
    delta_vel_m = delta_vel_m[~numpy.isnan(delta_vel_m)]
    delta_vel_over_tol = (abs(delta_vel_m)>self.binary_diff_tol)
    if (delta_vel_over_tol.any() == True):
      print " !!! warning: relative velocity difference after gravity_subset higher then", self.binary_diff_tol, "!!! "
      print "\t binary velocity vector differs by:"
      print "\t", delta_vel_m[delta_vel_over_tol]
      if diff_flag==0:
        diff_flag = 2
      elif diff_flag==1:
        diff_flag=3
    
    return diff_flag
    
  def evolve_model(self, end_time):
    self.particles = self.run_instances_in_serial(end_time)
    self.model_time = end_time
    
  def stop(self):
    for i,instance_i in enumerate(self.list_of_instances):
      instance_i.stop()
  
  def add_particles(self):
    for i,instance_i in enumerate(self.list_of_instances):
      instance_i.particles.add_particles(self.list_binary_and_stones_subsets[i])
  
  def commit_particles(self):
    for i,instance_i in enumerate(self.list_of_instances):
      instance_i.commit_particles()
  
  def initialize_code(self):
    """
    initialize n_subsets of code instances
    with parameters specified in options
    """
    list_of_instances = []
    for i in range(self.number_of_subsets):
      instance_i = self.gravity_code(self.converter, channel_type="sockets")
      instance_i.initialize_code()
      for parameter_name, parameter_value in self.options.iteritems():
        setattr(instance_i.parameters, parameter_name, parameter_value)
      instance_i.commit_parameters()
      list_of_instances.append(instance_i)
    return list_of_instances
  
  def kinetic_energy(self):
    """
    kinetic energy of the first subset (also taken as the result of evolve)
    (assuming the binary is composed from the only two particles with mass,
     this should be the same for all subsets)
    """
    return self.list_of_instances[0].kinetic_energy
  
  def potential_energy(self):
    """
    potential energy of the first subset (also taken as the result of evolve)
    (assuming the binary is composed from the only two particles with mass,
     this should be the same for all subsets)
    """
    return self.list_of_instances[0].potential_energy

    
  