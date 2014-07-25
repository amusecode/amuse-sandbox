from amuse.community.seba.interface import SeBa
from amuse.units import units, constants
from amuse.datamodel import Particles
from amuse.support.console import set_printing_strategy

from binary import *
from math import sqrt
import numpy as np

from amuse.community.seculartriple_TPS.interface import SecularTriple

from amuse.units import quantities
import matplotlib.pyplot as plt
import amuse.plot as aplt

#constants
time_step_factor = 0.01

# lowering this to 0.005 makes the code twice as slow
# 0.01 -> error in the semi-major axis of about 0.5%
maximum_wind_mass_loss_factor = 0.01 
error_dm = 0.05
maximum_radius_change_factor = 0.05 
error_dr = 0.2
minimum_time_step = 1.e-9 |units.Myr

REPORT_TRIPLE_EVOLUTION = False

class Triple:
    #-------
    #setup stellar system
    def __init__(self, inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            mutual_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node,
            metallicity,
            tend):
            
        if inner_semimajor_axis >= outer_semimajor_axis:
            print 'error input parameters, should be:'
            print 'inner_semimajor_axis < outer_semimajor_axis' 
            exit(-1)            
            
        stars = self.make_stars(inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis)
        bins = self.make_bins(stars, inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node)

        triples = Particles(1)
        triples[0].child1 = bins[0]
        triples[0].child2 = bins[1]
        triples[0].mass = triples[0].child1.mass + triples[0].child2.mass
        #when going to quadruples, mutual_inclination should be a characteristic of a bin 
        triples[0].mutual_inclination = mutual_inclination 
        triples[0].is_star = False
        triples[0].is_binary = True # True if there are 2 children
        triples[0].is_container = True # Upper level
        
        self.first_contact = True
        self.tend = tend #...
        self.time = 0.0|units.yr
        self.previous_time = 0.0|units.yr

        self.particles = triples
        self.setup_se_code(metallicity, stars)
        self.setup_secular_code(triples)

        self.update_previous_se_parameters()
        self.update_se_parameters() 
        self.update_se_wind_parameters()

    def make_stars(self, inner_primary_mass, inner_secondary_mass, outer_mass, inner_semimajor_axis, outer_semimajor_axis):
        stars = Particles(3)
        stars.is_star = True
        stars.is_binary = False 
        stars.is_container = False
        stars.is_donor = False

        stars[0].mass = inner_primary_mass
        stars[1].mass = inner_secondary_mass
        stars[2].mass = outer_mass

        # default now corotating
        corotating_angular_frequency_inner = corotating_spin_angular_frequency_binary(inner_semimajor_axis, stars[0].mass, stars[1].mass)
        corotating_angular_frequency_outer = corotating_spin_angular_frequency_binary(outer_semimajor_axis, stars[0].mass+stars[1].mass, stars[2].mass)

        stars[0].spin_angular_frequency = corotating_angular_frequency_inner 
        stars[1].spin_angular_frequency = corotating_angular_frequency_inner 
        stars[2].spin_angular_frequency = corotating_angular_frequency_outer

        return stars 
         
    def make_bins(self, stars, inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node):     
        bins = Particles(2)
        bins.is_star = False
        bins.is_binary = True
        bins.is_container = False
        bins.is_stable = True
        

        bins[0].child1 = stars[0]
        bins[0].child2 = stars[1]
        bins[0].child1.parent = bins[0]
        bins[0].child2.parent = bins[0]

        bins[0].semimajor_axis = inner_semimajor_axis
        bins[0].eccentricity = inner_eccentricity
        bins[0].argument_of_pericenter = inner_argument_of_pericenter
        bins[0].longitude_of_ascending_node = inner_longitude_of_ascending_node
        
        bins[0].mass_transfer_rate = 0.0 | units.MSun/units.yr
        bins[0].accretion_efficiency_mass_transfer = 1.0
        bins[0].accretion_efficiency_wind_child1_to_child2 = 0.0
        bins[0].accretion_efficiency_wind_child2_to_child1 = 0.0

        bins[1].child1 = stars[2]
        bins[1].child2 = bins[0]
        bins[1].child1.parent = bins[1]
        bins[1].child2.parent = bins[1]
        
        bins[1].semimajor_axis = outer_semimajor_axis
        bins[1].eccentricity = outer_eccentricity
        bins[1].argument_of_pericenter = outer_argument_of_pericenter                
        bins[1].longitude_of_ascending_node = outer_longitude_of_ascending_node
        
        bins[1].mass_transfer_rate = 0.0 | units.MSun/units.yr        
        bins[1].accretion_efficiency_mass_transfer = 1.0
        bins[1].accretion_efficiency_wind_child1_to_child2 = 0.0
        bins[1].accretion_efficiency_wind_child2_to_child1 = 0.0


        # binary evolutionary settings
        bins[0].specific_AM_loss_mass_transfer = 2.5 
        bins[1].specific_AM_loss_mass_transfer = 2.5

        
        for x in bins:
            x.mass = x.child1.mass + x.child2.mass

        return bins
    #-------
            
    #-------
    #setup community codes
    def setup_se_code(self, metallicity, stars):
        self.se_code = SeBa()
        self.se_code.parameters.metallicity = metallicity
        self.se_code.particles.add_particles(stars)
        self.channel_from_se = self.se_code.particles.new_channel_to(stars)
        self.channel_to_se = stars.new_channel_to(self.se_code.particles)
        self.channel_from_se.copy()
    
    def setup_secular_code(self, triples):
        self.secular_code = SecularTriple()
        self.secular_code.triples.add_particles(triples)
        self.secular_code.parameters.equations_of_motion_specification = 0
        self.secular_code.parameters.include_quadrupole_terms = True
        self.secular_code.parameters.include_octupole_terms = True        
        self.secular_code.parameters.include_inner_tidal_terms = False #need apsidal motion constant
        self.secular_code.parameters.include_outer_tidal_terms = False
        self.secular_code.parameters.include_inner_wind_terms = True
        self.secular_code.parameters.include_outer_wind_terms = True
        self.secular_code.parameters.include_inner_RLOF_terms = False
        self.secular_code.parameters.include_outer_RLOF_terms = False
        self.secular_code.parameters.include_1PN_inner_terms = False
        self.secular_code.parameters.include_1PN_outer_terms = False
        self.secular_code.parameters.include_1PN_inner_outer_terms = False ### warning: probably broken
        self.secular_code.parameters.include_25PN_inner_terms = False
        self.secular_code.parameters.include_25PN_outer_terms = False

        self.secular_code.parameters.check_for_dynamical_stability = True
        self.secular_code.parameters.check_for_inner_collision = True
        self.secular_code.parameters.check_for_outer_collision = False
        self.secular_code.parameters.check_for_inner_RLOF = False ### work in progress
        self.secular_code.parameters.check_for_outer_RLOF = False ### work in progress

        self.channel_from_secular = self.secular_code.triples.new_channel_to(triples)
        self.channel_to_secular = triples.new_channel_to(self.secular_code.triples)
    #-------

    #-------
    def update_previous_se_parameters(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.particles[0]

        self.previous_time = self.time

        if stellar_system.is_container:
            self.update_previous_se_parameters(stellar_system.child2)
            stellar_system.previous_mass = stellar_system.mass 
        elif stellar_system.is_star:
            stellar_system.previous_mass = stellar_system.mass      
            stellar_system.previous_radius = stellar_system.radius
        elif stellar_system.is_binary:
            self.update_previous_se_parameters(stellar_system.child1)        
            self.update_previous_se_parameters(stellar_system.child2)
            stellar_system.previous_mass = stellar_system.mass 
        else:
            print 'update_previous_se_parameters: structure stellar system unknown'        
            exit(-1)
            
#        self.particles[0].child1.child1.previous_mass = self.particles[0].child1.child1.mass 
#        self.particles[0].child1.child2.previous_mass = self.particles[0].child1.child2.mass 
#        self.particles[0].child2.child1.previous_mass = self.particles[0].child2.child1.mass 
#        self.particles[0].child2.child2.previous_mass = self.particles[0].child2.child2.mass 

#        self.particles[0].child1.child1.previous_radius = self.particles[0].child1.child1.radius 
#        self.particles[0].child1.child2.previous_radius = self.particles[0].child1.child2.radius 
#        self.particles[0].child2.child1.previous_radius = self.particles[0].child2.child1.radius 
    #-------

    #-------
    def update_se_wind_parameters(self):
        self.update_wind_mass_loss_rate()
        self.update_time_derivative_of_radius()
                    
    def update_wind_mass_loss_rate(self, stellar_system = None):
        #note: wind mass loss rate < 0
        if stellar_system == None:
            stellar_system = self.particles[0]

        if stellar_system.is_container:
            self.update_wind_mass_loss_rate(stellar_system.child2)
        elif stellar_system.is_star:
            star_in_se_code = stellar_system.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
            stellar_system.wind_mass_loss_rate = star_in_se_code.get_wind_mass_loss_rate() 
        elif stellar_system.is_binary:
            self.update_wind_mass_loss_rate(stellar_system.child1)        
            self.update_wind_mass_loss_rate(stellar_system.child2)
        else:
            print 'update_wind_mass_loss_rate: structure stellar system unknown'        
            exit(-1)
        
#        ch1ch1_in_se_code = self.particles[0].child1.child1.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
#        ch1ch2_in_se_code = self.particles[0].child1.child2.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
#        ch2ch1_in_se_code = self.particles[0].child2.child1.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
#            
#        self.particles[0].child1.child1.wind_mass_loss_rate = ch1ch1_in_se_code.get_wind_mass_loss_rate()
#        self.particles[0].child1.child2.wind_mass_loss_rate = ch1ch2_in_se_code.get_wind_mass_loss_rate()
#        self.particles[0].child2.child1.wind_mass_loss_rate = ch2ch1_in_se_code.get_wind_mass_loss_rate()

           
    def update_time_derivative_of_radius(self, stellar_system = None):
        #update time_derivative_of_radius for effect of wind on spin
        #radius change due to stellar evolution, not mass transfer
        if stellar_system == None:
            stellar_system = self.particles[0]
                
        time_step = self.time - self.previous_time

        if self.time == quantities.zero:
            #initialization
            self.particles[0].child1.child1.time_derivative_of_radius = 0.0 | units.RSun/units.yr
            self.particles[0].child1.child2.time_derivative_of_radius = 0.0 | units.RSun/units.yr
            self.particles[0].child2.child1.time_derivative_of_radius = 0.0 | units.RSun/units.yr
        else:    
            if stellar_system.is_container:
                self.update_time_derivative_of_radius(stellar_system.child2)
            elif stellar_system.is_star:
                stellar_system.time_derivative_of_radius = (stellar_system.radius - stellar_system.previous_radius)/time_step
            elif stellar_system.is_binary:
                self.update_time_derivative_of_radius(stellar_system.child1)        
                self.update_time_derivative_of_radius(stellar_system.child2)
            else:
                print 'update_time_derivative_of_radius: structure stellar system unknown'        
                exit(-1)


#        if time_step > 0|units.yr:
#            self.particles[0].child1.child1.time_derivative_of_radius = (self.particles[0].child1.child1.radius - self.particles[0].child1.child1.previous_radius)/time_step
#            self.particles[0].child1.child2.time_derivative_of_radius = (self.particles[0].child1.child2.radius - self.particles[0].child1.child2.previous_radius)/time_step
#            self.particles[0].child2.child1.time_derivative_of_radius = (self.particles[0].child2.child1.radius - self.particles[0].child2.child1.previous_radius)/time_step
#        else:
            #initialization
#            self.particles[0].child1.child1.time_derivative_of_radius = 0.0 | units.RSun/units.yr
#            self.particles[0].child1.child2.time_derivative_of_radius = 0.0 | units.RSun/units.yr
#            self.particles[0].child2.child1.time_derivative_of_radius = 0.0 | units.RSun/units.yr
    #-------

    #-------
    def update_se_parameters(self):
        self.update_binary_mass()
        self.update_stellar_parameters()
        
    def update_binary_mass(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.particles[0]

        if stellar_system.is_container:
            self.update_binary_mass(stellar_system.child2)
            stellar_system.mass = stellar_system.child1.mass + stellar_system.child2.mass
        elif stellar_system.is_star:
            return
        elif stellar_system.is_binary:
            self.update_binary_mass(stellar_system.child1)        
            self.update_binary_mass(stellar_system.child2)
            stellar_system.mass = stellar_system.child1.mass + stellar_system.child2.mass
        else:
            print 'update_binary_mass: structure stellar system unknown'        
            exit(-1)
      
#        if self.particles[0].child1.is_binary:
#            self.particles[0].child1.mass = self.particles[0].child1.child1.mass + self.particles[0].child1.child2.mass
#        if self.particles[0].child2.is_binary:
#            self.particles[0].child2.mass = self.particles[0].child2.child1.mass + self.particles[0].child2.child2.mass
#    

    def update_stellar_parameters(self, stellar_system = None):
        # for the convective envelope mass:
        # the prescription of Hurley, Pols & Tout 2000 is implemented in SeBa, however note that the prescription in BSE is different
        # for the convective envelope radius:
        # the prescription of Hurley, Tout & Pols 2002 is implemented in SeBa, however note that the prescription in BSE is different
        # both parameters don't need to be updated manually anymore
        
        
        if stellar_system == None:
            stellar_system = self.particles[0]

        if stellar_system.is_container:
            self.update_stellar_parameters(stellar_system.child2)
        elif stellar_system.is_star:
            star_in_se_code = stellar_system.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
            stellar_system.gyration_radius = star_in_se_code.get_gyration_radius_sq()**0.5     
            stellar_system.apsidal_motion_constant = 1. # WARNING
        elif stellar_system.is_binary:
            self.update_stellar_parameters(stellar_system.child1)        
            self.update_stellar_parameters(stellar_system.child2)
        else:
            print 'update_previous_se_parameters: structure stellar system unknown'        
            exit(-1)
        
#        self.particles[0].child1.child1.convective_envelope_mass = self.particles[0].child1.child1.convective_envelope_mass
#        self.particles[0].child1.child2.convective_envelope_mass = self.particles[0].child1.child2.convective_envelope_mass
#        self.particles[0].child2.child1.convective_envelope_mass = self.particles[0].child2.child1.convective_envelope_mass
#
#        self.particles[0].child1.child1.convective_envelope_radius = self.particles[0].child1.child1.convective_envelope_radius
#        self.particles[0].child1.child2.convective_envelope_radius = self.particles[0].child1.child2.convective_envelope_radius
#        self.particles[0].child2.child1.convective_envelope_radius = self.particles[0].child2.child1.convective_envelope_radius
#
#        self.particles[0].child1.child1.gyration_radius = self.se_code.particles[0].get_gyration_radius_sq()**0.5
#        self.particles[0].child1.child2.gyration_radius = self.se_code.particles[1].get_gyration_radius_sq()**0.5
#        self.particles[0].child2.child1.gyration_radius = self.se_code.particles[2].get_gyration_radius_sq()**0.5
    #-------

    #-------
    # useful functions general
    
    #whether or not a binary system consists of just two stars
    def is_double_star(self, stellar_system=None):
        if stellar_system == None:
            stellar_system = self.particles[0]    
    
        if stellar_system.is_binary and stellar_system.child1.is_star and stellar_system.child2.is_star:
            return True
        else:
            return False

    def is_triple(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.particles[0]

        # this structure is not completely nice, for quintiples also child1 should be checked
        if stellar_system.is_binary and stellar_system.child2.is_binary:
            if self.is_double_star(stellar_system.child2.child2) and stellar_system.child2.child1.is_star:
                return True
            elif self.is_double_star(stellar_system.child2.child1) and stellar_system.child2.child2.is_star:
                print 'is_triple structure possible?' 
                return True

        return False

    def has_donor(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.particles[0]
            
        if stellar_system.is_container:
            if self.has_donor(stellar_system.child2):
                return True
        elif stellar_system.is_star:
            if stellar_system.is_donor:
                return True
        elif stellar_system.is_binary:
            if self.has_donor(stellar_system.child1) or self.has_donor(stellar_system.child2):
                return True                        
        else:
            print 'has_donor: structure stellar system unknown'        
            exit(-1)
            
        return False            

    def kozai_timescale(self):
        if self.is_triple():
            P_in = orbital_period(self.particles[0].child2.child2) #period inner binary 
            P_out = orbital_period(self.particles[0].child2)#period outer binary 
            #Kozai 1962, Michaealy & Perets 2014        
#            return P_out**2 / P_in * (self.particles[0].child2.child2.mass / self.particles[0].child2.child1.mass)
            #Kinoshita & Nakai 1999, Hamers et al 2014
            alpha_kozai = 1.
            return alpha_kozai * P_out**2 / P_in * (self.particles[0].child2.mass / self.particles[0].child2.child1.mass) * (1-self.particles[0].child2.eccentricity**2)**1.5       
        else:
            print 'Kozai timescale needs triple system'
            return np.nan   
            
    def check_for_RLOF(self):
        if self.is_triple():

            Rl1, Rl2, Rl3 = self.secular_code.give_roche_radii(self.particles[0])
            
            if self.particles[0].child1.child1.radius >= Rl1:
                self.particles[0].child1.child1.is_donor = True
            else:                                                         
               self.particles[0].child1.child1.is_donor = False

            if self.particles[0].child1.child2.radius >= Rl2:
                self.particles[0].child1.child2.is_donor = True
            else:                                                         
               self.particles[0].child1.child2.is_donor = False

            if self.particles[0].child2.child1.radius >= Rl3:
                self.particles[0].child2.child1.is_donor = True
            else:                                                         
               self.particles[0].child2.child1.is_donor = False
        elif self.is_double_star():
            Rl1 = roche_radius(self, self.child1)
            Rl2 = roche_radius(self, self.child2)
            
            if self.particles[0].child1.radius >= Rl1:
                self.particles[0].child1.is_donor = True
            else:                                                         
               self.particles[0].child1.child1.is_donor = False

            if self.particles[0].child2.radius >= Rl2:
                self.particles[0].child2.is_donor = True
            else:                                                         
               self.particles[0].child1.child2.is_donor = False
        else:
            print 'print_stellar_sytem: structure stellar system unknown'        
            exit(-1)           
    #-------

    #-------
    # useful functions for printing
    def print_star(self, star):
        if star.is_star:
            print 'star:'
            print star.age, 
            print star.stellar_type, 
            print star.mass, 
            print star.radius, 
            print star.core_mass, 
            print star.core_radius,
            print star.envelope_mass,
            print star.convective_envelope_mass,
            print star.convective_envelope_radius,
            print star.CO_core_mass,
            print star.luminosity,
            print star.temperature,
            print star.wind_mass_loss_rate
            print star.spin_angular_frequency,
            print star.is_donor
            print '\t'             
        else:
            print 'print_star needs a star'        
            exit(-1)
    
    
    def print_binary(self, binary):
        if binary.is_binary:
            print binary.mass, 
            print binary.semimajor_axis, 
            print binary.eccentricity, 
            print binary.argument_of_pericenter, 
            print binary.longitude_of_ascending_node,
            print binary.mass_transfer_rate,
#            print binary.mass_transfer_timescale,
            print binary.accretion_efficiency_mass_transfer,
            print binary.accretion_efficiency_wind_child1_to_child2,
            print binary.accretion_efficiency_wind_child2_to_child1,
            print binary.specific_AM_loss_mass_transfer,
            print binary.is_stable
            print '\t'
        else:
            print 'print_binary needs a binary'        
            exit(-1)
    
    
    def print_stellar_system(self, stellar_system):
        if stellar_system.is_container:
            print '\t'
            print 'stellar system: '
            print stellar_system.mutual_inclination
            print '\t'
            self.print_stellar_system(stellar_system.child2)
        elif stellar_system.is_binary:
            print 'binary star: '
            self.print_binary(stellar_system)
            self.print_stellar_system(stellar_system.child1)
            self.print_stellar_system(stellar_system.child2)
        elif stellar_system.is_star:
            self.print_star(stellar_system)
        else:
            print 'print_stellar_sytem: structure stellar system unknown'        
            exit(-1)
        print '\t'            
    #-------
        
    #-------
    # determining time steps
    def determine_time_step_wind(self, stellar_system = None):
    #note: returned value can be inf when the wind_mass_loss_rate = 0
        if stellar_system == None:
            stellar_system = self.particles[0]
            
        if stellar_system.is_container:
            dt = self.determine_time_step_wind(stellar_system.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_wind_triple = ", dt
            return dt
        elif stellar_system.is_star:
            dt = np.inf |units.Myr
            if stellar_system.wind_mass_loss_rate * -1. > quantities.zero:
               dt = maximum_wind_mass_loss_factor*stellar_system.mass / stellar_system.wind_mass_loss_rate*-1
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_wind_star = ", dt
            return dt 
        elif stellar_system.is_binary:
            dt1 = self.determine_time_step_wind(stellar_system.child1)        
            dt2 = self.determine_time_step_wind(stellar_system.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_wind_binary = ", dt1, dt2
            return min(dt1, dt2) 
        else:
            print 'determine_time_step_wind: structure stellar system unknown'        
            exit(-1)
    
    def determine_time_step_radius_change(self, stellar_system = None):
    #note: returned value can be inf when the change in radius <= 0
        #radius is only necessary for tides
        if not self.secular_code.parameters.include_inner_tidal_terms and not self.secular_code.parameters.include_outer_tidal_terms:
            return np.inf |units.Myr
    
        if stellar_system == None:
            stellar_system = self.particles[0]
            
        if stellar_system.is_container:
            dt = self.determine_time_step_radius_change(stellar_system.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_radius_change_triple = ", dt
            return dt
        elif stellar_system.is_star:
            dt = np.inf |units.Myr
            if stellar_system.time_derivative_of_radius > quantities.zero:
               dt = maximum_radius_change_factor*stellar_system.radius / stellar_system.time_derivative_of_radius
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_radius_change_star = ", dt
            return dt 
        elif stellar_system.is_binary:
            dt1 = self.determine_time_step_radius_change(stellar_system.child1)        
            dt2 = self.determine_time_step_radius_change(stellar_system.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_radius_change_binary = ", dt1, dt2
            return min(dt1, dt2) 
        else:
            print 'determine_time_step_radius_change: structure stellar system unknown'        
            exit(-1)
    
 
 
    def determine_time_step_mt(self, stellar_system = None):
    #note: returned value can be inf when the mass_transfer_rate = 0 
        if stellar_system == None:
            stellar_system = self.particles[0]
            
        if stellar_system.is_container:
            dt = self.determine_time_step_mt(stellar_system.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_mt_triple = ", dt
            return dt
        elif stellar_system.is_star:
            dt = np.inf |units.Myr
            if stellar_system.is_donor:
                dt = abs(time_step_factor*stellar_system.mass/stellar_system.mass_transfer_rate)
                print 'mt >0?: ', stellar_system.mass_transfer_rate
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_mt_star = ", dt
            return dt 
        elif stellar_system.is_binary:
            if stellar_system.child1.is_donor and stellar_system.child2.is_donor:
                print 'determine_time_step_mt: contact system'        
                exit(-1)

            dt1 = self.determine_time_step_mt(stellar_system.child1)        
            dt2 = self.determine_time_step_mt(stellar_system.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_mt_binary = ", dt1, dt2
            return min(dt1, dt2) 
        else:
            print 'determine_time_step_mt: structure stellar system unknown'        
            exit(-1)
            
    
    def determine_time_step(self):         
        if REPORT_TRIPLE_EVOLUTION:
            print "Dt = ", self.se_code.particles.time_step, self.tend/100.0
    
        #during unstable mass transfer or other instantaneous interactions
        if not self.particles[0].child2.is_stable:
            #no step in time
            return
        if not self.particles[0].child1.is_stable:
            # no step in time
            return
            
        #maximum time_step            
        time_step = self.tend - self.time

        # time_step of stellar evolution
        time_step =  min(time_step, min(self.se_code.particles.time_step))    
#        if self.particles[0].child2.child1.is_star:
#            time_step = min(time_step, self.particles[0].child2.child1.time_step)
#        if self.particles[0].child1.child1.is_star:
#            time_step = min(time_step, self.particles[0].child1.child1.time_step)
#        if self.particles[0].child1.child2.is_star:
#            time_step = min(time_step, self.particles[0].child1.child2.time_step)
                    
        # small time_step during heavy wind mass losses
        time_step = min(time_step, self.determine_time_step_wind())

#        if -1*self.particles[0].child1.child1.wind_mass_loss_rate * time_step / self.particles[0].child1.child1.mass > maximum_wind_mass_loss_factor:
#            time_step = min(time_step, maximum_wind_mass_loss_factor*self.particles[0].child1.child1.mass / self.particles[0].child1.child1.wind_mass_loss_rate*-1)
#        if -1*self.particles[0].child1.child2.wind_mass_loss_rate * time_step / self.particles[0].child1.child2.mass > maximum_wind_mass_loss_factor:
#            time_step = min(time_step, maximum_wind_mass_loss_factor*self.particles[0].child1.child2.mass / self.particles[0].child1.child2.wind_mass_loss_rate*-1)
#        if -1*self.particles[0].child2.child1.wind_mass_loss_rate * time_step / self.particles[0].child2.child1.mass > maximum_wind_mass_loss_factor:
#            time_step = min(time_step, maximum_wind_mass_loss_factor*self.particles[0].child2.child1.mass / self.particles[0].child2.child1.wind_mass_loss_rate*-1)
                         
        # small time_step during phases of fast growth (note: not yet during fast shrinkage)
        time_step = min(time_step, self.determine_time_step_radius_change())           
            
        #during stable mass transfer     
        time_step = min(time_step, self.determine_time_step_mt())

#        if (self.particles[0].child1.child1.is_donor):
#            time_step = min(time_step, abs(time_step_factor*self.particles[0].child1.child1.mass/self.particles[0].child1.mass_transfer_rate))
#        if (self.particles[0].child1.child2.is_donor):
#            time_step = min(time_step, abs(time_step_factor*self.particles[0].child1.child2.mass/self.particles[0].child1.mass_transfer_rate))
#        if (self.particles[0].child2.child1.is_donor):
#            time_step = min(time_step, abs(time_step_factor*self.particles[0].child2.child1.mass/self.particles[0].child2.mass_transfer_rate))
            
            
        ### for testing/plotting purposes only ###
#        time_step = min(time_step, self.tend/100.0)
    
            
            
        if time_step < minimum_time_step:
            print 'error small time_step'
            exit(-1)    
        time_step = max(time_step, minimum_time_step)                
        self.time += time_step
    #-------
        
    #-------
    #evolution
    def resolve_interaction_in_stellar_system(self):
    
        if REPORT_TRIPLE_EVOLUTION:
            print '\ninner binary'
        if self.is_double_star(self.particles[0].child1):
            resolve_binary_interaction(self.particles[0].child1, self)
        elif self.particles[0].child1.is_star:
    #        'e.g. if system merged'
            print 'do nothing'
        else:
            print 'resolve triple interaction: type of inner system unknown'
            exit(-1)                    
    
    
        if REPORT_TRIPLE_EVOLUTION:
            print '\nouter binary'
        elif self.particles[0].child2.is_binary:
            resolve_binary_interaction(self.particles[0].child2, self)
        else:
            print 'resolve triple interaction: type of outer system unknown'
            exit(-1)  
            
            
    def determine_mass_transfer_timescale(self):
    
        if self.particles[0].child2.child1.is_donor:
            if self.particles[0].child1.child1.is_donor or self.particles[0].child1.child2.is_donor:
                print 'RLOF in inner and outer binary'
                exit(0)
    
        if REPORT_TRIPLE_EVOLUTION:
            print '\ninner binary'
        if self.is_double_star(self.particles[0].child1):
            mass_transfer_stability(self.particles[0].child1)
        elif self.particles[0].child1.is_star:
            print 'do nothing'
        else:
            print 'determine_mass_transfer_timescale: type of inner system unknown'
            exit(-1)                    
    
        if REPORT_TRIPLE_EVOLUTION:
            print '\nouter binary'
        if self.particles[0].child2.is_binary:
            mass_transfer_stability(self.particles[0].child2)
        else:
            print 'determine_mass_transfer_timescale: type of outer system unknown'
            exit(-1)           


    def evolve_triple(self):
        
        # for plotting data
        times_array = quantities.AdaptingVectorQuantity() 
        a_in_array = quantities.AdaptingVectorQuantity()
        e_in_array = []
        g_in_array = []
        o_in_array = []
        a_out_array = quantities.AdaptingVectorQuantity()
        e_out_array = []
        g_out_array = []
        o_out_array = []
        i_mutual_array = []
        m1_array = quantities.AdaptingVectorQuantity()
        m2_array = quantities.AdaptingVectorQuantity()
        m3_array = quantities.AdaptingVectorQuantity()
    
        times_array.append(self.time)
        e_in_array.append(self.particles[0].child1.eccentricity)
        a_in_array.append(self.particles[0].child1.semimajor_axis)
        g_in_array.append(self.particles[0].child1.argument_of_pericenter) 
        o_in_array.append(self.particles[0].child1.longitude_of_ascending_node)               
        e_out_array.append(self.particles[0].child2.eccentricity)
        a_out_array.append(self.particles[0].child2.semimajor_axis)
        g_out_array.append(self.particles[0].child2.argument_of_pericenter) 
        o_out_array.append(self.particles[0].child2.longitude_of_ascending_node)               
        i_mutual_array.append(self.particles[0].mutual_inclination)        
        m1_array.append(self.particles[0].child1.child1.mass)
        m2_array.append(self.particles[0].child1.child2.mass)
        m3_array.append(self.particles[0].child2.child1.mass)

        print 'kozai timescale:', self.kozai_timescale()      
        while self.time<self.tend:
            self.update_previous_se_parameters()
#            self.determine_mass_transfer_timescale()
            self.determine_time_step()  
    
            #do stellar evolution 
            self.se_code.evolve_model(self.time)
            self.channel_from_se.copy()
            self.update_se_wind_parameters()
            self.update_se_parameters()        
            safety_check_time_step(self)
        
            # do secular evolution
            self.channel_to_secular.copy()   
            if self.is_triple() == True:
                self.secular_code.evolve_model(self.time)
            else:# e.g. binaries
                print 'Secular code disabled'
                exit(-1)
    
            if self.secular_code.triples[0].dynamical_instability == True:
                print "Dynamical instability at time/Myr = ",self.time.value_in(units.Myr)
                exit(0)
            if self.secular_code.triples[0].inner_collision == True:
                print "Inner collision at time/Myr = ",self.time.value_in(units.Myr)
                exit(0)
            if self.secular_code.triples[0].outer_collision == True:
                print "Outer collision at time/Myr = ",self.time.value_in(units.Myr)
                exit(0)
            self.channel_from_secular.copy()     
    
            #What if only 1 star left...
            self.check_for_RLOF()

            # when the secular code finds that mass transfer starts, go back in time
            if self.has_donor() and self.first_contact:
                    if REPORT_TRIPLE_EVOLUTION:
                        print 'Times:', self.previous_time, self.time, self.secular_code.model_time
                    print 'RLOF not yet'
#                    self.determine_mass_transfer_timescale()
                    break

                        
#                    self.time = self.secular_code.model_time 
#                    
#                    # mass should not be transferred just yet -> check if this works ok
#                    # do not overwrite this parameter in the secular code    
#                    self.particles[0].child1.child1.is_donor = False
#                    self.particles[0].child1.child2.is_donor = False
#                    self.particles[0].child2.child1.is_donor = False
#                    
#                    self.first_contact = False
#                    
    
           
            self.resolve_interaction_in_stellar_system()        
            #should also do safety check time_step here
            
            
            
            
            
            
            
            
            # for plotting data
            times_array.append(self.time)
            e_in_array.append(self.particles[0].child1.eccentricity)
            a_in_array.append(self.particles[0].child1.semimajor_axis)
            g_in_array.append(self.particles[0].child1.argument_of_pericenter) 
            o_in_array.append(self.particles[0].child1.longitude_of_ascending_node)               
            e_out_array.append(self.particles[0].child2.eccentricity)
            a_out_array.append(self.particles[0].child2.semimajor_axis)
            g_out_array.append(self.particles[0].child2.argument_of_pericenter) 
            o_out_array.append(self.particles[0].child2.longitude_of_ascending_node)               
            i_mutual_array.append(self.particles[0].mutual_inclination)        
            m1_array.append(self.particles[0].child1.child1.mass)
            m2_array.append(self.particles[0].child1.child2.mass)
            m3_array.append(self.particles[0].child2.child1.mass)
            
            
        # for plotting data
        e_in_array = np.array(e_in_array)
        g_in_array = np.array(g_in_array)
        o_in_array = np.array(o_in_array)
        e_out_array = np.array(e_out_array)
        g_out_array = np.array(g_out_array)
        o_out_array = np.array(o_out_array)
        i_mutual_array = np.array(i_mutual_array)

        self.plot_data = plot_data_container()
        self.plot_data.times_array = times_array
        self.plot_data.a_in_array = a_in_array
        self.plot_data.e_in_array = e_in_array
        self.plot_data.g_in_array = g_in_array
        self.plot_data.o_in_array = o_in_array        
        self.plot_data.a_out_array = a_out_array
        self.plot_data.e_out_array = e_out_array
        self.plot_data.g_out_array = g_out_array
        self.plot_data.o_out_array = o_out_array        
        self.plot_data.i_mutual_array = i_mutual_array
        self.plot_data.m1_array = m1_array
        self.plot_data.m2_array = m2_array
        self.plot_data.m3_array = m3_array
    #-------



def safety_check_time_step(triple):
        dm_1 = (triple.particles[0].child1.child1.previous_mass - triple.particles[0].child1.child1.mass)/triple.particles[0].child1.child1.mass
        dm_2 = (triple.particles[0].child1.child2.previous_mass - triple.particles[0].child1.child2.mass)/triple.particles[0].child1.child2.mass
        dm_3 = (triple.particles[0].child2.child1.previous_mass - triple.particles[0].child2.child1.mass)/triple.particles[0].child2.child1.mass
                
        if REPORT_TRIPLE_EVOLUTION:
            print 'wind mass loss rates:', 
            print triple.particles[0].child1.child1.wind_mass_loss_rate,
            print triple.particles[0].child1.child2.wind_mass_loss_rate,
            print triple.particles[0].child2.child1.wind_mass_loss_rate
            print 'relative wind mass losses:',
            print dm_1, dm_2, dm_3
        if (dm_1 > error_dm) or (dm_2 > error_dm) or (dm_3 > error_dm):
            print 'Change in mass in a single time_step larger then', error_dm
            print dm_1, dm_2, dm_3
            print triple.particles[0].child1.child1.stellar_type
            print triple.particles[0].child1.child2.stellar_type
            print triple.particles[0].child2.child1.stellar_type
#            print 'WARNING'
            exit(-1)


        if triple.secular_code.parameters.include_inner_tidal_terms or triple.secular_code.parameters.include_outer_tidal_terms:
            dr_1 = (triple.particles[0].child1.child1.radius - triple.particles[0].child1.child1.previous_radius)/triple.particles[0].child1.child1.radius
            dr_2 = (triple.particles[0].child1.child2.radius - triple.particles[0].child1.child2.previous_radius)/triple.particles[0].child1.child2.radius
            dr_3 = (triple.particles[0].child2.child1.radius - triple.particles[0].child2.child1.previous_radius)/triple.particles[0].child2.child1.radius
    
    
            if REPORT_TRIPLE_EVOLUTION:    
                print 'change in radius over time:', 
                print triple.particles[0].child1.child1.time_derivative_of_radius,
                print triple.particles[0].child1.child2.time_derivative_of_radius,
                print triple.particles[0].child2.child1.time_derivative_of_radius
                print 'relative change in radius:',
                print dr_1, dr_2, dr_3
            if (dr_1 > error_dr) or (dr_2 > error_dr) or (dr_3 > error_dr):
                print 'Change in radius in a single time_step larger then', error_dr
                print dr_1, dr_2, dr_3
                print triple.particles[0].child1.child1.stellar_type
                print triple.particles[0].child1.child2.stellar_type
                print triple.particles[0].child2.child1.stellar_type
    #            print 'WARNING'
                exit(-1)
            

    

class plot_data_container():
    def __init__(self):
        return

def plot_function(triple):
    ### plots to test secular code ###

    times_array_Myr = triple.plot_data.times_array.value_in(units.Myr)
    t_max_Myr = max(times_array_Myr)
    a_in_array_AU = triple.plot_data.a_in_array.value_in(units.AU)
    g_in_array = triple.plot_data.g_in_array
    e_in_array = triple.plot_data.e_in_array
    i_mutual_array = triple.plot_data.i_mutual_array
    o_in_array = triple.plot_data.o_in_array
    a_out_array_AU = triple.plot_data.a_out_array.value_in(units.AU)
    g_out_array = triple.plot_data.g_out_array
    e_out_array = triple.plot_data.e_out_array
    o_out_array = triple.plot_data.o_out_array
    m1_array = triple.plot_data.m1_array.value_in(units.MSun)
    m2_array = triple.plot_data.m2_array.value_in(units.MSun)
    m3_array = triple.plot_data.m3_array.value_in(units.MSun)
    
    figure = plt.figure(figsize=(10,13))
    N_subplots = 4

    plot_e_in = figure.add_subplot(N_subplots,1,1)
    plot_i_mutual = figure.add_subplot(N_subplots,1,2)
    plot_e_in_g_in = figure.add_subplot(N_subplots,1,3)
    plot_a_in = figure.add_subplot(N_subplots,1,4)


    plot_e_in.plot(times_array_Myr,e_in_array)
    plot_e_in.set_xlim(0,t_max_Myr)
    plot_e_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_e_in.set_ylabel('$e_\mathrm{in}$')

    plot_i_mutual.plot(times_array_Myr,i_mutual_array*180.0/np.pi)
    plot_i_mutual.set_xlim(0,t_max_Myr)
    plot_i_mutual.set_ylim(0.9*min(i_mutual_array*180.0/np.pi),1.1*max(i_mutual_array*180.0/np.pi))
    plot_i_mutual.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_mutual.set_ylabel('$i_\mathrm{mutual} ({}^\circ)$')

    plot_e_in_g_in.plot(np.cos(g_in_array),e_in_array)
    plot_e_in_g_in.set_xlabel('$\cos(g_\mathrm{in})$')
    plot_e_in_g_in.set_ylabel('$e_\mathrm{in}$')

    plot_a_in.plot(times_array_Myr,a_in_array_AU)
    plot_a_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_in.set_ylabel('$a_\mathrm{in}$')
    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)

    generic_plot_name = '_M'+str(m1_array[0]) + '_m'+str(m2_array[0]) +'_n'+str(m3_array[0]) + '_a'+str(a_in_array_AU[0]) + '_A'+str(a_out_array_AU[0]) + '_e'+str(e_in_array[0]) + '_E'+str(e_out_array[0]) + '_i'+str(i_mutual_array[0]/np.pi*180.0) + '_g'+str(g_in_array[0]) + '_G'+str(g_out_array[0]) + '_o'+str(o_in_array[0]) + '_O'+str(o_out_array[0]) + '_t'+str(t_max_Myr)
    plt.savefig('plots/orbit/TPS_inner_orbit'+generic_plot_name+'.pdf')
    plt.show()



    Mtot = m1_array+m2_array    
    plt.plot(times_array_Myr,a_in_array_AU)
    plt.plot(times_array_Myr,a_in_array_AU, '.')
    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot)
    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot, '.')
    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*Mtot[:-1]/Mtot[1:])
    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*Mtot[:-1]/Mtot[1:], '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$a_\mathrm{in}$')
    plt.savefig('plots/orbit/semi_inner'+generic_plot_name+'.pdf')
    plt.show()
    
    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot/a_in_array_AU)
    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot/a_in_array_AU, '.')
    plt.ylabel('$relative error a_\mathrm{in}$')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.savefig('plots/orbit/semi_inner_rel'+generic_plot_name+'.pdf')
    plt.show()

    dm = (m1_array[1:] - m1_array[:-1] )
    dt = (times_array_Myr[1:] - times_array_Myr[:-1])
    dmdt = (m1_array[1:] - m1_array[:-1] )/(times_array_Myr[1:] - times_array_Myr[:-1])


    plt.plot(dmdt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:])
    plt.plot(dmdt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:], '.')
    plt.show()

    plt.plot(dm, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:])
    plt.plot(dm, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:], '.')
    plt.show()

    plt.plot(dt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:])
    plt.plot(dt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:], '.')
    plt.show()

    
    
    
    
    
#    plt.plot(times_array_Myr, g_in_array*180.0/np.pi%360)
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$g_\mathrm{in}$')
#    plt.show()
#    plt.plot(times_array_Myr, np.cos(g_in_array))
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$\cos(g_\mathrm{in})$')
#    plt.show()
#
#
#    plt.plot(times_array_Myr, o_in_array*180.0/np.pi%360)
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$o_\mathrm{in}$')
#    plt.show()
#    plt.plot(times_array_Myr, np.cos(o_in_array))
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$\cos(o_\mathrm{in})$')
#    plt.show()


    #outer binary
    figure = plt.figure(figsize=(10,13))
    N_subplots = 4

    plot_e_out = figure.add_subplot(N_subplots,1,1)
    plot_i_mutual2 = figure.add_subplot(N_subplots,1,2)
    plot_e_out_g_out = figure.add_subplot(N_subplots,1,3)
    plot_a_out = figure.add_subplot(N_subplots,1,4)

#    times_array_Myr = triple.plot_data.times_array.value_in(units.Myr)
#    t_max_Myr = max(times_array_Myr)
#    i_mutual_array = triple.plot_data.i_mutual_array

    plot_e_out.plot(times_array_Myr,e_out_array)
    plot_e_out.set_xlim(0,t_max_Myr)
    plot_e_out.set_xlabel('$t/\mathrm{Myr}$')
    plot_e_out.set_ylabel('$e_\mathrm{out}$')

    plot_i_mutual2.plot(times_array_Myr,i_mutual_array*180.0/np.pi)
    plot_i_mutual2.set_xlim(0,t_max_Myr)
    plot_i_mutual2.set_ylim(0.9*min(i_mutual_array*180.0/np.pi),1.1*max(i_mutual_array*180.0/np.pi))
    plot_i_mutual2.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_mutual2.set_ylabel('$i_\mathrm{mutual} ({}^\circ)$')

    plot_e_out_g_out.plot(np.cos(g_out_array),e_out_array)
    plot_e_out_g_out.set_xlabel('$\cos(g_\mathrm{out})$')
    plot_e_out_g_out.set_ylabel('$e_\mathrm{out}$')

    plot_a_out.plot(times_array_Myr,a_out_array_AU)
    plot_a_out.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_out.set_ylabel('$a_\mathrm{out}$')

    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)
    plt.savefig('plots/orbit/TPS_outer_orbit'+generic_plot_name+'.pdf')
    plt.show()


    Mtott = m1_array+m2_array+m3_array    
    plt.plot(times_array_Myr,a_out_array_AU)
    plt.plot(times_array_Myr,a_out_array_AU, '.')
    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott)
    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott, '.')
    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*Mtott[:-1]/Mtott[1:])
    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*Mtott[:-1]/Mtott[1:], '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$a_\mathrm{out}$')
    plt.savefig('plots/orbit/semi_outer'+generic_plot_name+'.pdf')
    plt.show()

    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott/a_out_array_AU)
    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott/a_out_array_AU, '.')
    plt.ylabel('$relative error a_\mathrm{out}$')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.savefig('plots/orbit/semi_outer_rel'+generic_plot_name+'.pdf')
    plt.show()

    dm = (m3_array[1:] - m3_array[:-1] )
    dmdt = (m3_array[1:] - m3_array[:-1] )/(times_array_Myr[1:] - times_array_Myr[:-1])

    plt.plot(dmdt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:])
    plt.plot(dmdt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:], '.')
    plt.show()

    plt.plot(dm, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:])
    plt.plot(dm, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:], '.')
    plt.show()

    plt.plot(dt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:])
    plt.plot(dt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:], '.')
    plt.show()



#    plt.plot(np.cos(g_out_array), e_in_array)
#    plt.xlabel('$\cos(g_\mathrm{out})$')
#    plt.ylabel('$e_\mathrm{in}$')
#    plt.show()
#
#    plt.plot(np.cos(g_in_array), np.cos(g_out_array))
#    plt.xlabel('$\cos(g_\mathrm{in})$')
#    plt.ylabel('$\cos(g_\mathrm{out})$')
#    plt.show()
#
#    plt.plot(times_array_Myr, g_out_array*180.0/np.pi%360)
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$g_\mathrm{out}$')
#    plt.show()
#    plt.plot(times_array_Myr, np.cos(g_out_array))
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$\cos(g_\mathrm{out})$')
#    plt.show()
#
#    plt.plot(times_array_Myr, o_out_array*180.0/np.pi%360)
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$o_\mathrm{out}$')
#    plt.show()
#    plt.plot(times_array_Myr, np.cos(o_out_array))
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$\cos(o_\mathrm{out})$')
#    plt.show()
#
#    aplt.plot(times_array_Myr, m1_array)
#    aplt.plot(times_array_Myr, m1_array, '.')
#    aplt.plot(times_array_Myr, m2_array)
#    aplt.plot(times_array_Myr, m2_array, '.')
#    aplt.plot(times_array_Myr, m3_array)
#    aplt.plot(times_array_Myr, m3_array, '.')
#    plt.xlim(0, t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$M/\mathrm{MSun}$')
#    plt.show()
#    
    
#-----
#for running triple.py from other routines
def main(inner_primary_mass= 1.3|units.MSun, inner_secondary_mass= 0.5|units.MSun, outer_mass= 0.5|units.MSun,
            inner_semimajor_axis= 1.0 |units.AU, outer_semimajor_axis= 100.0 |units.AU,
            inner_eccentricity= 0.1, outer_eccentricity= 0.5,
            mutual_inclination= 80.0*np.pi/180.0,
            inner_argument_of_pericenter= 0.1, outer_argument_of_pericenter= 0.5,
            inner_longitude_of_ascending_node= 0.0, outer_longitude_of_ascending_node= 0.0,
            metallicity= 0.02|units.none,
            tend= 5.0 |units.Myr):

    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

            


    inner_eccentricity = float(inner_eccentricity)
    outer_eccentricity = float(outer_eccentricity)
    mutual_inclination = float(mutual_inclination)
    inner_argument_of_pericenter = float(inner_argument_of_pericenter)
    outer_argument_of_pericenter = float(outer_argument_of_pericenter)
    inner_longitude_of_ascending_node = float(inner_longitude_of_ascending_node)
    outer_longitude_of_ascending_node = float(outer_longitude_of_ascending_node)


    triple = Triple(inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            mutual_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node,
            metallicity,
            tend)

    triple.evolve_triple()
#    plot_function(triple)
    triple.print_stellar_system(triple.particles[0])
    return triple
#-----

#-----
#for running triple.py from the commandline
def parse_arguments():
    from amuse.units.optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-a", unit=units.AU,
                      dest="inner_semimajor_axis", type="float", 
                      default = 1.0 |units.AU,
                      help="inner semi major axis [%default]")
    parser.add_option("-A", unit=units.AU,
                      dest="outer_semimajor_axis", type="float", 
                      default = 100.0 |units.AU,
                      help="outer semi major axis [%default]")
    parser.add_option("-e",
                      dest="inner_eccentricity", type="float", default = 0.1,
                      help="inner eccentricity [%default]")
    parser.add_option("-E",
                      dest="outer_eccentricity", type="float", default = 0.5,
                      help="outer eccentricity [%default]")
    parser.add_option("-i",
                      dest="mutual_inclination", type="float", default = 80.0*np.pi/180.0,
                      help="mutual inclination [rad] [%default]")
    parser.add_option("-g",
                      dest="inner_argument_of_pericenter", type="float", default = 0.1,
                      help="inner argument of pericenter [rad] [%default]")
    parser.add_option("-G",
                      dest="outer_argument_of_pericenter", type="float", default = 0.5,
                      help="outer argument of pericenter [rad] [%default]")
    parser.add_option("-o",
                      dest="inner_longitude_of_ascending_node", type="float", default = 0.0,
                      help="inner longitude of ascending node [rad] [%default]")
    parser.add_option("-O",
                      dest="outer_longitude_of_ascending_node", type="float", default = 0.0,
                      help="outer longitude of ascending node [rad] [%default]")
    parser.add_option("-M", unit=units.MSun, 
                      dest="inner_primary_mass", type="float", default = 1.3|units.MSun,
                      help="inner primary mass [%default]")
    parser.add_option("-m",  unit=units.MSun, 
                      dest="inner_secondary_mass", type="float", default = 0.5|units.MSun,
                      help="inner secondary mass [%default]")
    parser.add_option("-l",  unit=units.MSun, 
                      dest="outer_mass", type="float", default = 0.5|units.MSun,
                      help="outer mass [%default]")
    parser.add_option("-t", "-T", unit=units.Myr, 
                      dest="tend", type="float", default = 5.0 |units.Myr,
                      help="end time [%default] %unit")
    parser.add_option("-z", unit=units.none, 
                      dest="metallicity", type="float", default = 0.02|units.none,
                      help="metallicity [%default] %unit")

    options, args = parser.parse_args()
    return options.__dict__





if __name__ == '__main__':
    options = parse_arguments()

    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")


    triple = Triple(**options)
    triple.evolve_triple()
    plot_function(triple)
