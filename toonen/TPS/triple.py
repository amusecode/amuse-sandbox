## Triple:      Triple evolution
##              computes the evolution of a given triple
##              given any initial conditions (M, m, l, A, a, E, e, i, G, g, O, o).

## Output:      in the form of the following files:
##              

## Options:   -M    mass of the inner primary [Msun]
##            -m    mass of the inner secondary [MSun]
##            -l    mass of the outer star [Msun]

##            -A    inner semi-major axis []
##            -a    outer semi-major axis []

##            -E    inner eccentricity [1.]
##            -e    outer eccentricity [0.]

##            -i    relative inclination between orbits []

##            -G    inner argument of pericenter []
##            -g    outer argument of pericenter []

##            -O    inner longitude of ascending nodes []
##            -o    outer longitude of ascending nodes []

##            -T or -t   binary end time. [13500 Myr]
##            -z         metallicity of stars  [0.02 Solar] 


from amuse.community.seba.interface import SeBa
from binary import *
import sys 
sys.path.insert(1, '../../hamers/TPS/code')
from seculartriple_TPS.interface import SecularTriple

from amuse.units import units, constants
from amuse.datamodel import Particles
from amuse.support.console import set_printing_strategy
from amuse.io import write_set_to_file
from amuse.units import quantities
import amuse.plot as aplt

import matplotlib.pyplot as plt
from math import sqrt
import numpy as np


REPORT_TRIPLE_EVOLUTION = False 

#stopping conditions
#for the moment these are global variables, and not interface parameters, as we cannot simulate them
stop_at_merger = True # as implementation is missing, e.g. secular code not adjusted
stop_at_disintegrated = True # as implementation is missing, e.g. secular code not adjusted
stop_at_triple_mass_transfer = True # as implementation is missing 
stop_at_collision = True # as implementation is missing
stop_at_dynamical_instability = True # this should always be true!


file_name = "triple.hdf"
file_type = "hdf5"
#file_name = "triple.txt"
#file_type = "txt"

#constants
time_step_factor_stable_mt = 0.01 #1% mass loss during mass transfer
# lowering this to 0.005 makes the code twice as slow
# 0.01 -> error in the semi-major axis of about 0.5%
maximum_wind_mass_loss_factor = 0.01 
error_dm = 0.05
maximum_radius_change_factor = 0.005
error_dr = 0.05
minimum_time_step = 1.e-9 |units.Myr
min_mass = 0.08 |units.MSun # for stars
max_mass = 100 |units.MSun
stellar_types_SN_remnants = [13,14]|units.stellar_type # remnant types created through a supernova
minimum_outer_eccentricity = 1.e-3

class Triple_Class:
    #-------
    #setup stellar system
    def __init__(self, inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            relative_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node,
            metallicity, tend, number):      
            
        outer_eccentricity = self.test_initial_parameters(inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis, inner_eccentricity, outer_eccentricity,
            relative_inclination, inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node)   
        if inner_primary_mass < inner_secondary_mass:
            spare = inner_primary_mass
            inner_primary_mass = inner_secondary_mass
            inner_secondary_mass = spare     
                        
        stars = self.make_stars(inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis)
        bins = self.make_bins(stars, inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node)
        
        self.first_contact = True # not used at the moment, how to reset ?
        self.instantaneous_evolution = False # no secular evolution        
        self.tend = tend #...
        self.time = 0.0|units.yr
        self.previous_time = 0.0|units.yr

        self.triple = bins[1]
        self.triple.relative_inclination = relative_inclination 
        self.triple.is_star = False
        self.triple.dynamically_stable = True 
        self.triple.number = number 
        
        self.setup_stellar_code(metallicity, stars)
        self.setup_secular_code(self.triple.as_set())

        self.update_previous_stellar_parameters()
        self.update_stellar_parameters() 
        self.update_stellar_wind_parameters()


    def make_stars(self, inner_primary_mass, inner_secondary_mass, outer_mass, inner_semimajor_axis, outer_semimajor_axis):
        stars = Particles(3)
        stars.is_star = True
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
        bins.is_stable = True
        bins.part_dt_mt = 1.
        bins.bin_type = bin_type['unknown'] #Unknown

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

        return bins
    #-------
            
    #-------
    #setup community codes
    def setup_stellar_code(self, metallicity, stars):
        self.stellar_code = SeBa()

        #stopping conditions:
#        print self.stellar_code.stopping_conditions
#        print self.stellar_code.stopping_conditions.supernova_detection.is_supported()
#        print self.stellar_code.stopping_conditions.supernova_detection.is_enabled()
#        print self.stellar_code.stopping_conditions.supernova_detection.is_set()
#        print self.stellar_code.stopping_conditions.supernova_detection.enable()
#        print self.stellar_code.stopping_conditions.supernova_detection.is_enabled()
#        print self.stellar_code.stopping_conditions.supernova_detection.is_set()

        self.stellar_code.parameters.metallicity = metallicity
        self.stellar_code.particles.add_particles(stars)
        self.channel_from_stellar = self.stellar_code.particles.new_channel_to(stars)
        self.channel_to_stellar = stars.new_channel_to(self.stellar_code.particles)
        self.channel_from_stellar.copy()

      
    def setup_secular_code(self, triple_set):
        self.secular_code = SecularTriple()
        self.secular_code.triples.add_particles(triple_set)
        self.secular_code.parameters.equations_of_motion_specification = 0
        self.secular_code.parameters.include_quadrupole_terms = True
        self.secular_code.parameters.include_octupole_terms = True        
        self.secular_code.parameters.include_inner_tidal_terms = False
        self.secular_code.parameters.include_outer_tidal_terms = False
        self.secular_code.parameters.include_inner_wind_terms = True
        self.secular_code.parameters.include_outer_wind_terms = True
        self.secular_code.parameters.include_inner_RLOF_terms = True
        self.secular_code.parameters.include_outer_RLOF_terms = True
        self.secular_code.parameters.include_1PN_inner_terms = False
        self.secular_code.parameters.include_1PN_outer_terms = False
        self.secular_code.parameters.include_1PN_inner_outer_terms = False ### warning: probably broken
        self.secular_code.parameters.include_25PN_inner_terms = False
        self.secular_code.parameters.include_25PN_outer_terms = False

        self.secular_code.parameters.check_for_dynamical_stability = True
        self.secular_code.parameters.check_for_inner_collision = True
        self.secular_code.parameters.check_for_outer_collision = True
        self.secular_code.parameters.check_for_inner_RLOF = False ### work in progress
        self.secular_code.parameters.check_for_outer_RLOF = False ### work in progress

        self.channel_from_secular = self.secular_code.triples.new_channel_to(triple_set)
        self.channel_to_secular = triple_set.new_channel_to(self.secular_code.triples)
    #-------

    #-------
    def test_initial_parameters(self, inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            relative_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node):
            

        if max(inner_primary_mass, outer_mass) > max_mass:  
#        if (min(inner_secondary_mass, outer_mass) < min_mass) or (max(inner_primary_mass, outer_mass) > max_mass):  
            print inner_primary_mass, inner_secondary_mass, outer_mass
            print 'should be within:', min_mass, '-', max_mass
            print 'error: masses not in allowed range'
            exit(1)
                        
        if inner_semimajor_axis >= outer_semimajor_axis:
            print 'error input parameters, should be:'
            print 'inner_semimajor_axis < outer_semimajor_axis' 
            exit(1)            
        if (inner_semimajor_axis < 0.|units.RSun):
            print 'error: inner separation not in allowed range'
            exit(1)
        if (outer_semimajor_axis < 0.|units.RSun):
            print 'error: outer separation not in allowed range'
            exit(1)
    
        if (inner_eccentricity < 0.) or (inner_eccentricity > 1.):
            print 'error: inner eccentricity not in allowed range'
            exit(1)
        if (outer_eccentricity < 0.) or (outer_eccentricity > 1.):
            print 'error: outer eccentricity not in allowed range'
            exit(1)
        if (outer_eccentricity < minimum_outer_eccentricity):
            outer_eccentricity = minimum_outer_eccentricity
    
        if (relative_inclination < 0.) or (relative_inclination > 2.*np.pi):
            print 'error: relative inclination not in allowed range'
            exit(1)
    
        if (inner_argument_of_pericenter < 0.) or (inner_argument_of_pericenter > 2*np.pi):
            print 'error: inner argument of pericenter not in allowed range'
            exit(1)
        if (outer_argument_of_pericenter < 0.) or (outer_argument_of_pericenter > 2*np.pi):
            print 'error: outer argument of pericenter not in allowed range'
            exit(1)
    
        if (inner_longitude_of_ascending_node < 0.) or (inner_longitude_of_ascending_node > 2*np.pi):
            print 'error: inner longitude of ascending node not in allowed range'
            exit(1)
        if (outer_longitude_of_ascending_node < 0.) or (outer_longitude_of_ascending_node > 2*np.pi):
            print 'error: outer longitude of ascending node not in allowed range'
            exit(1)
            
        return outer_eccentricity            
    #-------

    #-------
    def update_previous_stellar_parameters(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        self.previous_time = self.time
        if stellar_system.is_star:
            stellar_system.previous_mass = self.get_mass(stellar_system)      
            stellar_system.previous_radius = stellar_system.radius
            stellar_system.previous_stellar_type = stellar_system.stellar_type
        else:
            self.update_previous_stellar_parameters(stellar_system.child1)        
            self.update_previous_stellar_parameters(stellar_system.child2)
            stellar_system.previous_mass = self.get_mass(stellar_system) 
    #-------

    #-------
    def update_stellar_wind_parameters(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple
            
        self.update_wind_mass_loss_rate(stellar_system)
        self.update_time_derivative_of_radius(stellar_system)
                    
    def update_wind_mass_loss_rate(self, stellar_system = None):
        #note: wind mass loss rate < 0
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
#            if stellar_system.mass > 0.15|units.MSun:
#                stellar_system.radius = 0.88824945031|units.RSun
#            elif stellar_system.mass < 0.09|units.MSun:
#                stellar_system.radius = 0.14498995417 |units.RSun
#            else:
#                stellar_system.radius = 0.13045599412  |units.RSun   
#            

            star_in_stellar_code = stellar_system.as_set().get_intersecting_subset_in(self.stellar_code.particles)[0]
            stellar_system.wind_mass_loss_rate = star_in_stellar_code.get_wind_mass_loss_rate() 
        else:
            self.update_wind_mass_loss_rate(stellar_system.child1)        
            self.update_wind_mass_loss_rate(stellar_system.child2)
        
           
    def update_time_derivative_of_radius(self, stellar_system = None):
        #update time_derivative_of_radius for effect of wind on spin
        #radius change due to stellar evolution, not mass transfer
        if stellar_system == None:
            stellar_system = self.triple
                
        time_step = self.time - self.previous_time

        if self.time == quantities.zero:
            #initialization
            self.triple.child2.child1.time_derivative_of_radius = 0.0 | units.RSun/units.yr
            self.triple.child2.child2.time_derivative_of_radius = 0.0 | units.RSun/units.yr
            self.triple.child1.time_derivative_of_radius = 0.0 | units.RSun/units.yr
        else:     
            if stellar_system.is_star:
                stellar_system.time_derivative_of_radius = (stellar_system.radius - stellar_system.previous_radius)/time_step
            else:
                self.update_time_derivative_of_radius(stellar_system.child1)        
                self.update_time_derivative_of_radius(stellar_system.child2)
    #-------

    #-------
    def update_stellar_parameters(self, stellar_system = None):
        # for the convective envelope mass:
        # the prescription of Hurley, Pols & Tout 2000 is implemented in SeBa, however note that the prescription in BSE is different
        # for the convective envelope radius:
        # the prescription of Hurley, Tout & Pols 2002 is implemented in SeBa, however note that the prescription in BSE is different
        # both parameters don't need to be updated manually anymore
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            star_in_stellar_code = stellar_system.as_set().get_intersecting_subset_in(self.stellar_code.particles)[0]
            stellar_system.gyration_radius = star_in_stellar_code.get_gyration_radius_sq()**0.5     
            stellar_system.apsidal_motion_constant = self.apsidal_motion_constant(stellar_system) 
            if stellar_system.convective_envelope_radius < 0|units.RSun:
                print 'convective_envelope_radius < 0'
                exit(1)
            if stellar_system.convective_envelope_radius == 0|units.RSun:
                stellar_system.convective_envelope_mass = 1.e-10 |units.MSun    
                stellar_system.convective_envelope_radius = 1.e-10 |units.RSun    
        else:
            self.update_stellar_parameters(stellar_system.child1)        
            self.update_stellar_parameters(stellar_system.child2)
    #-------

    #-------
    # useful functions general
    
    #whether or not a stellar system consists of just two stars
    def is_binary(self, stellar_system=None):
        if stellar_system == None:
            stellar_system = self.triple    
    
        if not stellar_system.is_star and stellar_system.child1.is_star and stellar_system.child2.is_star:
            return True
        else:
            return False

    def is_triple(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if not stellar_system.is_star:
            if stellar_system.child1.is_star and self.is_binary(stellar_system.child2):
                return True
            elif stellar_system.child2.is_star and self.is_binary(stellar_system.child2):
                return True

        return False

    def has_donor(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            if stellar_system.is_donor:
                return True
        else:
            if self.has_donor(stellar_system.child1) or self.has_donor(stellar_system.child2):
                return True                        
            
        return False            


    def has_contact_system(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            return False
        elif self.is_binary(stellar_system):
            if stellar_system.child1.is_donor and stellar_system.child2.is_donor:
                return True
        else:
            if self.has_contact_system(stellar_system.child1):
                return True
            if self.has_contact_system(stellar_system.child2):
                return True
            
        return False            

# if a merger is currently taking place, not if a merger has happened in the past
    def has_merger(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            return False
        else:
            if self.has_merger(stellar_system.child1):
                return True
            if self.has_merger(stellar_system.child2):
                return True
            if stellar_system.bin_type == bin_type['merger']:  
                return True    

        return False            

# if a disruption is currently taking place, not if a disruption has happened in the past
    def has_disintegrated(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            return False
        else:
            if self.has_disintegrated(stellar_system.child1):
                return True
            if self.has_disintegrated(stellar_system.child2):
                return True
            if stellar_system.bin_type == bin_type['disintegrated']:  
                return True    
            
        return False            

# if a mass transfer in the outer binary of the triple is currently taking place, not if a mass transfer has happened in the past
    def has_triple_mass_transfer(self, stellar_system = None): 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            return False
        elif self.is_binary(stellar_system):
            return False
        else:
            if self.has_triple_mass_transfer(stellar_system.child1):
                return True
            if self.has_triple_mass_transfer(stellar_system.child2):
                return True
            if stellar_system.bin_type != bin_type['unknown'] and stellar_system.bin_type != bin_type['detached']:  
                return True    
            
        return False            


    def has_stellar_type_changed(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if stellar_system.stellar_type != stellar_system.previous_stellar_type:
                return True
        else:
            if self.has_stellar_type_changed(stellar_system.child1) or self.has_stellar_type_changed(stellar_system.child2):
                return True                        
            
        return False            
    


    def is_system_stable(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            return True
        elif self.is_binary(stellar_system):
            return stellar_system.is_stable
        else:
            if stellar_system.is_stable and self.is_system_stable(stellar_system.child1) and self.is_system_stable(stellar_system.child2):
                return True                        
            
        return False                    

    def get_mass(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            return stellar_system.mass
        else:
            M1 = self.get_mass(stellar_system.child1)        
            M2 = self.get_mass(stellar_system.child2)
            return M1 + M2 
    #-------

    #-------
    # useful functions general
            
    def orbital_period(self, bs):
        Porb = 2*np.pi * np.sqrt(bs.semimajor_axis**3/constants.G / self.get_mass(bs))
        return Porb

    def orbital_angular_momentum(self, bs):
        M = self.get_mass(bs.child1)
        m = self.get_mass(bs.child2)
        a = bs.semimajor_axis
        e = bs.eccentricity
        J = M*m * np.sqrt(constants.G*a*(1-e**2)/(M+m))
    
        if REPORT_BINARY_EVOLUTION:
            print 'Jorb:', M, m, a, e, J
    
        return J
    
    def stellar_angular_momentum(self, ss):
        if ss.is_star:
            moment_of_inertia = ss.gyration_radius**2 * ss.mass * ss.radius**2
            Jstar = moment_of_inertia * ss.spin_angular_frequency
            return Jstar            
        else:
            print 'stellar_angular_momentum: structure stellar system unknown'        
            exit(2)
            
    def apsidal_motion_constant(self, star):
        if star.stellar_type in [13]|units.stellar_type: #ns
            #based on Brooke & Olle 1955, for n=1 polytrope
            return 0.260
    
        elif star.stellar_type in [1,7,10,11,12]|units.stellar_type:#ms, he-ms, wd
            #based on Brooke & Olle 1955, for n=3 polytrope
            return 0.0144            

        elif star.stellar_type in [0,2,3,4,5,6,8,9,17]|units.stellar_type:#low-mass ms, hg, gb, cheb, agb, he-g, pre-ms
            #based on Brooke & Olle 1955, for n=3 polytrope
#            return 0.143 
            #based on Claret & Gimenez 1992, 96, 225 the value should be smaller, try:
            return 0.05
        else:
            print 'apsidal motion constant: stellar_type unknown'
#            print 'possibly black hole'
            print star.stellar_type
            exit(2)
            

    def kozai_timescale(self):
        if self.is_triple():
           alpha_kozai = 1.
           if self.triple.child1.is_star:
                star = self.triple.child1
                bin = self.triple.child2
           else: 
                star = self.triple.child2
                bin = self.triple.child1


           P_in = self.orbital_period(bin) #period inner binary 
           P_out = self.orbital_period(self.triple)#period outer binary 
           return alpha_kozai * P_out**2 / P_in * (self.get_mass(self.triple) / self.get_mass(star)) * (1-self.triple.eccentricity**2)**1.5       

        else:
            print 'Kozai timescale needs triple system'
            return np.nan   
            
    def check_for_RLOF(self):
        if self.triple.is_star:
            return
        elif self.is_binary():
            Rl1 = roche_radius(self, self.child1)
            Rl2 = roche_radius(self, self.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print 'Roche lobe radii:', Rl1, Rl2
                print 'Stellar radii:', self.triple.child1.radius, self.triple.child2.radius
            
            if self.triple.child1.radius >= Rl1:
                self.triple.child1.is_donor = True
            else:                                                         
               self.triple.child1.child1.is_donor = False

            if self.triple.child2.radius >= Rl2:
                self.triple.child2.is_donor = True
            else:                                                         
               self.triple.child1.child2.is_donor = False
 
        elif self.is_triple():
            if self.triple.child1.is_star:
                star = self.triple.child1
                bin = self.triple.child2
            else:
                star = self.triple.child2
                bin = self.triple.child1

            #assumping secular code always returns inner binary first
            Rl1, Rl2, Rl3 = self.secular_code.give_roche_radii(self.triple)
            if REPORT_TRIPLE_EVOLUTION:
                print 'Roche lobe radii:', Rl1, Rl2, Rl3
                print 'Stellar radii:', bin.child1.radius, bin.child2.radius, star.radius

            bin.child1.is_donor = False
            bin.child2.is_donor = False
            star.is_donor = False
            
            if bin.child1.radius >= Rl1:
                bin.child1.is_donor = True
            if bin.child2.radius >= Rl2:
                bin.child2.is_donor = True
            if star.radius >= Rl3:
                star.is_donor = True
                
            if star.is_donor and (bin.child1.is_donor or bin.child2.is_donor):
                print 'RLOF in inner and outer binary'
                print Rl1, bin.child1.radius, Rl2, bin.child2.radius
                print RL3, star.radius
                exit(1)                   
                
        else:
            print 'check_for_RLOF: structure stellar system unknown'        
            exit(2)    
                     
#            
#    def determine_partial_timestep_stable_mass_transfer(self, stellar_system = None):
#        if stellar_system == None:
#            stellar_system = self.triple
#
#        if stellar_system.is_star:
#            return np.inf |units.Myr 
#        else:
#            dt1 = self.determine_partial_timestep_stable_mass_transfer(stellar_system.child1)        
#            dt2 = self.determine_partial_timestep_stable_mass_transfer(stellar_system.child2)
#            dt =  stellar_system.part_dt_mt
#            return min(dt, min(dt1, dt2))
                   
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
            print star.wind_mass_loss_rate,
            print star.spin_angular_frequency,
            print star.is_donor
            print '\t'             
        else:
            print 'print_star needs a star'        
            exit(2)
    
    
    def print_binary(self, binary):
        if not binary.is_star:
            print self.get_mass(binary), 
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
            exit(2)
    
    
    def print_stellar_system(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple
            print stellar_system.number,
            print stellar_system.relative_inclination

        if stellar_system.is_star:
            self.print_star(stellar_system)
        else:
            print 'binary star: '
            self.print_binary(stellar_system)
            self.print_stellar_system(stellar_system.child1)
            self.print_stellar_system(stellar_system.child2)
        print '\t'            
    #-------
    #-------
    #don't change this unless you know what you're doing
    def remove_parents(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        parents = []
        if stellar_system.is_star:
            try:
                p = stellar_system.parent
                stellar_system.parent = 0
                return p
            except AttributeError: #when there is no parent
                return 0
        else:
            parents.append(self.remove_parents(stellar_system.child1))
            parents.append(self.remove_parents(stellar_system.child2))

            p = stellar_system.parent
#            except AttributeError: #when there is no parent=
            if p != None:
                stellar_system.parent = 0
                parents.append(p)                                    

            return parents
            
    #don't change this unless you know what you're doing
    def set_parents(self, parents, stellar_system=None):            
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
                stellar_system.parent = parents
        else:
            self.set_parents(parents[0], stellar_system.child1)
            self.set_parents(parents[1], stellar_system.child2)
            if len(parents) == 3:
                stellar_system.parent = parents[2]
            elif len(parents) != 2:
                print 'set_parents: structure stellar system unknown'        
                exit(2)
 
    def save_snapshot(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if file_type == 'txt':
            parents = self.remove_parents()
            write_set_to_file(self.triple.as_set(), file_name, file_type) 
            self.set_parents(parents)
        else:
            write_set_to_file(self.triple.as_set(), file_name, file_type, version='2.0') 

    #some minor parameters are missing:
#        self.first_contact = True 
#        self.instantaneous_evolution = False 
#        self.tend = tend #...
#        self.time = 0.0|units.yr
#        self.previous_time = 0.0|units.yr
    #-------
        
    #-------
    # determining time steps
    def determine_time_step_wind(self, stellar_system = None):
    #note: returned value can be inf when the wind_mass_loss_rate = 0
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            dt = np.inf |units.Myr
            if stellar_system.wind_mass_loss_rate * -1. > quantities.zero:
               dt = maximum_wind_mass_loss_factor*stellar_system.mass / stellar_system.wind_mass_loss_rate*-1
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_wind_star = ", dt
            return dt 
        else:
            dt1 = self.determine_time_step_wind(stellar_system.child1)        
            dt2 = self.determine_time_step_wind(stellar_system.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_wind_binary = ", dt1, dt2
            return min(dt1, dt2) 
    
    def determine_time_step_radius_change(self, stellar_system = None):
    #note: returned value can be inf when the change in radius <= 0
        #radius is only necessary for tides
        if not self.secular_code.parameters.include_inner_tidal_terms and not self.secular_code.parameters.include_outer_tidal_terms:
            return np.inf |units.Myr
    
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            dt = np.inf |units.Myr
            if stellar_system.time_derivative_of_radius > quantities.zero:
               dt = maximum_radius_change_factor*stellar_system.radius / stellar_system.time_derivative_of_radius
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_radius_change_star = ", dt
            return dt 
        else:
            dt1 = self.determine_time_step_radius_change(stellar_system.child1)        
            dt2 = self.determine_time_step_radius_change(stellar_system.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_radius_change_binary = ", dt1, dt2
            return min(dt1, dt2) 
    
 
 
    def determine_time_step_stable_mt(self, stellar_system = None):
    #note: returned value can be inf when the mass_transfer_rate = 0 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            dt = np.inf |units.Myr
            if stellar_system.is_donor:
                dt = abs(time_step_factor_stable_mt*stellar_system.mass/stellar_system.parent.mass_transfer_rate)
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_mt_star = ", dt
            return dt
        else:
            if self.is_binary(stellar_system) and stellar_system.child1.is_donor and stellar_system.child2.is_donor:
                #should have been taken care of in determine_time_step()
                print 'determine_time_step_stable_mt: contact system'        
                exit(1)

            dt1 = self.determine_time_step_stable_mt(stellar_system.child1)        
            dt2 = self.determine_time_step_stable_mt(stellar_system.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print "Dt_mt_binary = ", dt1, dt2
            return min(dt1, dt2) 
            
    
    def determine_time_step(self):         
        if REPORT_TRIPLE_EVOLUTION:
            print "Dt = ", self.stellar_code.particles.time_step, self.tend/100.0
 
#        during unstable mass transfer, contact_system or other instantaneous interactions: no step in time
        if self.has_donor() and (not self.is_system_stable() or self.has_contact_system()):
            return            
                       
        #maximum time_step            
        time_step_max = self.tend - self.time
        
        # time_step of stellar evolution
        time_step_stellar_code = self.stellar_code.particles.time_step
                    
        # small time_step during heavy wind mass losses
        time_step_wind =  self.determine_time_step_wind()

        # small time_step during phases of fast growth (note: not yet during fast shrinkage)
        time_step_radius_change = self.determine_time_step_radius_change()             
        
        time_step = min(time_step_radius_change, min(time_step_wind, min( min(time_step_stellar_code), time_step_max)))    
#        print time_step_max, time_step_stellar_code, time_step_wind, time_step_radius_change, time_step#, time_step2  
          
        #during stable mass transfer     
        if self.has_donor():
            time_step = min(time_step, self.determine_time_step_stable_mt())

        ### for testing/plotting purposes only ###
#        time_step = min(time_step, self.tend/100.0)
#        time_step = min(time_step, min(time_step_stellar_code)/2.)
            
        if self.time == quantities.zero:
            #initialization (e.g. time_derivative_of_radius)
            #initially triple.child2 is the inner binary
            P_in = self.orbital_period(self.triple.child2) #period inner binary 
            time_step = min(P_in/10., time_step)

        if time_step < minimum_time_step:
            print 'error small time_step'
            print time_step_max, time_step_stellar_code, time_step_wind, time_step_radius_change, time_step
            exit(1)    
#        time_step = max(time_step, minimum_time_step)                
        self.time += time_step
    #-------

    #-------
    def adjust_system_after_supernova_kick(self):
        SN_star_in_stellar_code = self.stellar_code.stopping_conditions.supernova_detection.particles(0)
        print SN_star_in_stellar_code

        SN_star = SN_star_in_stellar_code.as_set().get_intersecting_subset_in(self.triple)[0] # this will probably not work :-)

        system = SN_star
        try:
            parent = system.parent             
        except AttributeError:            
            #SN_star is a single star
            #adjust velocity of SN_star
            # this should only be done for single stars, so for an exception in the first iteration
            return  
              
              
        while True:
            try:    
                system = system.parent
                self.adjust_binary_after_supernova_kick(system)  
            except AttributeError:
                #when there is no parent
                break


    def adjust_binary_after_supernova_kick(self, system):
        if self.is_binary(system):
            print 'adjust double star after supernova kick'
            #adjust orbit, separation and eccentricity
                       
            if system.eccentricity < 0:
                print 'e<0'
                exit(0)            
            elif system.eccentricity >= 1:
                print 'System becomes unbound'
                exit(0)
            else: 
                pericenter = system.semimajor_axis*(1-system.eccentricity)
                if pericenter < system.child1.radius + system.child2.radius:
                    print 'Collision'
                    exit(0)
            
            #adjust systematic velocity

        else: #binary    
            print 'adjust binary after supernova kick'
            #adjust orbit, separation and eccentricity, inclination
            #adjust systematic velocity

            if system.eccentricity < 0:
                print 'e<0'
                exit(0)            
            elif system.eccentricity >= 1:
                print 'System becomes unbound'
                exit(0)
            else: 
                pericenter = system.semimajor_axis*(1-system.eccentricity)
                if system.child1.is_star and not system.child2.is_star:
                    if pericenter < system.child1.radius + system.child2.semimajor_axis:
                        print 'Unstable?'
                        exit(0)
                elif not system.child1.is_star and system.child2.is_star:
                    if pericenter < system.child1.semimajor_axis + system.child2.radius:
                        print 'Unstable?'
                        exit(0)
                else:
                    print 'adjust_binary_after_supernova_kick: type of system unknown'
                    exit(2)
                
                

     
        
    #-------

    #-------
    #evolution
    def resolve_stellar_interaction(self, stellar_system = None):
    # the most inner binary is calculated first, and then move outwards

        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if REPORT_TRIPLE_EVOLUTION:
                print 'single stellar evolution'

            print 'for now no single stellar evolution - exiting program'
            exit(2)
                
            return
        elif self.is_binary(stellar_system):
            if REPORT_TRIPLE_EVOLUTION:
                print '\n resolve stellar interaction: binary - double star'            
            resolve_binary_interaction(stellar_system, self)
        else:
            if stellar_system.child2.is_star: #child1 is a binary
                self.resolve_stellar_interaction(stellar_system.child1)        
            elif stellar_system.child1.is_star: #child2 is a binary
                self.resolve_stellar_interaction(stellar_system.child2)        
            else:
                print 'resolve_stellar_interaction: structure stellar system unknown'        
                print 'both children are binaries'
                exit(2)
            
            if REPORT_TRIPLE_EVOLUTION:
                print '\n resolve stellar interaction: binary'
            resolve_binary_interaction(stellar_system, self)            
                            


#        if REPORT_TRIPLE_EVOLUTION:
#            print '\ninner binary - child1'
#        if self.triple.child1.is_star:
#            print 'do nothing'
#        elif self.is_double_star(self.triple.child1):
#            resolve_binary_interaction(self.triple.child1, self)
#        else:
#            print 'resolve triple interaction: type of inner system unknown'
#            exit(2)                    
        
#        if REPORT_TRIPLE_EVOLUTION:
#            print '\nouter binary - child2'
#        if self.triple.is_binary:
            #assumption that system has 3 stars or less -> only child1 and child2 exist
#            if self.triple.child2.is_binary:
#                resolve_binary_interaction(self.triple.child2, self)            
#            else:
#                print 'resolve triple interaction: type of outer system unknown'
#                exit(2)                    
                                    
    def determine_mass_transfer_timescale(self, stellar_system = None):

        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            if REPORT_TRIPLE_EVOLUTION:
                print 'single stellar evolution'
            return
        elif self.is_binary(stellar_system):
            if REPORT_TRIPLE_EVOLUTION:
                print '\n determine_mass_transfer_timescale: binary - double star'            
            mass_transfer_stability(stellar_system, self)
            if REPORT_TRIPLE_EVOLUTION:
                print 'mt rate double star:', stellar_system.mass_transfer_rate
        else:
            self.determine_mass_transfer_timescale(stellar_system.child1)        
            self.determine_mass_transfer_timescale(stellar_system.child2)        

            if REPORT_TRIPLE_EVOLUTION:
                print '\n determine_mass_transfer_timescale: binary'
            mass_transfer_stability(stellar_system, self)            
            if REPORT_TRIPLE_EVOLUTION:
                print 'mt rate binary:', stellar_system.mass_transfer_rate
                                
#   
#        if self.triple.is_binary and self.is_double_star(self.triple.child1) and self.triple.child2.is_binary:
#            #assumption child1 from self.triple.child2 is the star, and child2 is self.triple.child1
#            if self.triple.child2.child1.is_donor and (self.triple.child1.child1.is_donor or self.triple.child1.child2.is_donor):
#                print 'RLOF in inner and outer binary'
                
                
    
#        if REPORT_TRIPLE_EVOLUTION:
#            print '\ninner binary - child1'
#        if self.triple.child1.is_star:
#            print 'do nothing'
#        elif self.is_double_star(self.triple.child1):
#            mass_transfer_stability(self.triple.child1, self)
#            if REPORT_TRIPLE_EVOLUTION:
#                print 'mt rate double star:', self.triple.child1.mass_transfer_rate
#        else:
#            print 'resolve triple interaction: type of inner system unknown'
#            exit(2)                    
#        
#        if REPORT_TRIPLE_EVOLUTION:
#            print '\nouter binary - child2'
#        if self.triple.is_binary:
#            #assumption that system has 3 stars or less -> only child1 and child2 exist
#            if self.triple.child2.is_binary:
#                mass_transfer_stability(self.triple.child2, self)
#                if REPORT_TRIPLE_EVOLUTION:
#                    print 'mt rate binary:', self.triple.child2.mass_transfer_rate
#            else:
#                print 'resolve triple interaction: type of outer system unknown'
#                exit(2)                    
#    
    
    def safety_check_time_step(self, stellar_system = None):

        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            dm = (stellar_system.previous_mass - stellar_system.mass)/stellar_system.mass
            if REPORT_TRIPLE_EVOLUTION:
                print 'wind mass loss rate:', stellar_system.wind_mass_loss_rate,
                print 'relative wind mass losses:', dm

            if (dm > error_dm) and (not self.has_stellar_type_changed(stellar_system) or stellar_system.stellar_type not in stellar_types_SN_remnants):
                print 'Change in mass in a single time_step larger then', error_dm
                print dm, stellar_system.mass, stellar_system.stellar_type
                exit(1)

            if self.secular_code.parameters.include_inner_tidal_terms or self.secular_code.parameters.include_outer_tidal_terms:
                dr = (stellar_system.radius - stellar_system.previous_radius)/stellar_system.radius
                
                print 'radius:', stellar_system.radius,
                print 'change in radius over time:', stellar_system.time_derivative_of_radius,
                print 'relative change in radius:', dr
                print 'time_derivative_of_radius:', stellar_system.time_derivative_of_radius    
        
                if REPORT_TRIPLE_EVOLUTION:    
                    print 'change in radius over time:',  stellar_system.time_derivative_of_radius,
                    print 'relative change in radius:', dr
                if (dr > error_dr):
                    print 'Change in radius in a single time_step larger then', error_dr
                    print dr, stellar_system.mass, stellar_system.stellar_type
                    exit(1)

        else:
            self.safety_check_time_step(stellar_system.child1)        
            self.safety_check_time_step(stellar_system.child2)
            print 'eccentricity:', stellar_system.eccentricity
            
    
#-------------------------    

    def evolve_model(self):
        
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
        i_relative_array = []
        m1_array = quantities.AdaptingVectorQuantity()
        m2_array = quantities.AdaptingVectorQuantity()
        m3_array = quantities.AdaptingVectorQuantity()
    
        times_array.append(self.time)
        e_in_array.append(self.triple.child2.eccentricity)
        a_in_array.append(self.triple.child2.semimajor_axis)
        g_in_array.append(self.triple.child2.argument_of_pericenter) 
        o_in_array.append(self.triple.child2.longitude_of_ascending_node)               
        e_out_array.append(self.triple.eccentricity)
        a_out_array.append(self.triple.semimajor_axis)
        g_out_array.append(self.triple.argument_of_pericenter) 
        o_out_array.append(self.triple.longitude_of_ascending_node)               
        i_relative_array.append(self.triple.relative_inclination)        
        m1_array.append(self.triple.child2.child1.mass)
        m2_array.append(self.triple.child2.child2.mass)
        m3_array.append(self.triple.child1.mass)

        print 'kozai timescale:', self.kozai_timescale(), self.tend    
        self.determine_mass_transfer_timescale()
        self.save_snapshot()        
        while self.time<self.tend:
            if self.has_stellar_type_changed():
                self.save_snapshot()        
        
            self.update_previous_stellar_parameters()
            self.check_for_RLOF()   
            self.determine_time_step()  
            print '\n\ntime:', self.time, self.has_donor()            
    
            #do stellar evolution 
            self.stellar_code.evolve_model(self.time)
            self.channel_from_stellar.copy()
            
            self.update_stellar_wind_parameters()
            self.update_stellar_parameters()        
#            if  self.stellar_code.stopping_condition.supernova_detection.is_set():
#                print 'supernova detected'
#                print self.stellar_code.stopping_condition.supernova_detection.particles(0)
#                print self.stellar_code.stopping_condition.supernova_detection.particles(1)
#                self.adjust_system_after_supernova_kick()
#                exit(0)                   
#
#                print 'how do I restart the secular code?'                
#                self.secular_code.model_time = self.time
#                continue # the while loop, skip resolve_stellar_interaction and secular evolution
            
            #do stellar interaction
            self.determine_mass_transfer_timescale()
            self.resolve_stellar_interaction()
#            if  self.stellar_code.stopping_condition.supernova_detection.is_set():
#                print 'supernova detected'
#                print self.stellar_code.stopping_condition.supernova_detection.particles(0)
#                print self.stellar_code.stopping_condition.supernova_detection.particles(1)
#                self.adjust_system_after_supernova_kick()
#                exit(0)                   
#
#                print 'how do I restart the secular code?'                
#                self.secular_code.model_time = self.time
#                continue # the while loop, skip resolve_stellar_interaction and secular evolution

            if stop_at_merger and self.has_merger():
                print 'stop at merger in triple'
                break    
            if stop_at_disintegrated and self.has_disintegrated():
                print 'stopping at disintegration of triple'
                break
            if stop_at_triple_mass_transfer and self.has_triple_mass_transfer():
                print 'stopping at mass transfer in outer binary of triple'
                #it's possible that there is mass transfer in the inner and outer binary
                print self.triple.child2.bin_type
                print self.triple.bin_type
                print self.has_triple_mass_transfer()
                break                                    
                  

            # do secular evolution
            self.channel_to_secular.copy()   
            if self.instantaneous_evolution == False: 
                self.safety_check_time_step()
                if not self.is_triple:# e.g. binaries
                    print 'Secular code disabled'
                    exit(1)

                # if mass transfer should only take place for a fraction of the timestep
                # e.g. at the end of mass transfer when the envelope is thin
                if self.triple.child2.part_dt_mt < 1: # inner binary, see function determine_partial_time_step_stable_mass_transfer
                    full_dt = self.time - self.previous_time
                    self.secular_code.evolve_model(self.previous_time + full_dt * self.triple.child2.part_dt_mt)

                    if self.has_donor():
                        print 'After partial timestep the system should be detached...'
                        exit(1)
                    # not necessary because secular code reset is_donor and therefore the mass transfer rate is not used if the system is detached    
#                    self.triple.child2.mass_transfer_rate =  0.0 | units.MSun/units.yr 
                        
                    self.secular_code.evolve_model(self.time)
                    self.triple.child2.part_dt_mt = 1.                    
                    
                else:
                    self.secular_code.evolve_model(self.time)
                    
                    
                    
                if stop_at_dynamical_instability and self.secular_code.triples[0].dynamical_instability == True:
                    self.triple.dynamically_stable = False    
                    print "Dynamical instability at time/Myr = ",self.time.value_in(units.Myr)
                    self.save_snapshot()        
                    break
                if stop_at_collision and self.secular_code.triples[0].inner_collision == True:
                    self.triple.child2.bin_type = bin_type['collision']
                    print "Inner collision at time/Myr = ",self.time.value_in(units.Myr)
                    self.save_snapshot()        
                    break
                if stop_at_collision and self.secular_code.triples[0].outer_collision == True:
                    self.triple.bin_type = bin_type['collision']
                    print "Outer collision at time/Myr = ",self.time.value_in(units.Myr)
                    self.save_snapshot()        
                    break
                   
                self.channel_from_secular.copy()     
            else:        
                self.secular_code.model_time = self.time
                self.instantaneous_evolution = False
                
                 
            #should also do safety check time_step here -> only make sure that mass loss from stable mass transfer is not too large -> determine_time_step_mt
            

            
            # for plotting data
            times_array.append(self.time)
            e_in_array.append(self.triple.child2.eccentricity)
            a_in_array.append(self.triple.child2.semimajor_axis)
            g_in_array.append(self.triple.child2.argument_of_pericenter) 
            o_in_array.append(self.triple.child2.longitude_of_ascending_node)               
            e_out_array.append(self.triple.eccentricity)
            a_out_array.append(self.triple.semimajor_axis)
            g_out_array.append(self.triple.argument_of_pericenter) 
            o_out_array.append(self.triple.longitude_of_ascending_node)               
            i_relative_array.append(self.triple.relative_inclination)        
            m1_array.append(self.triple.child2.child1.mass)
            m2_array.append(self.triple.child2.child2.mass)
            m3_array.append(self.triple.child1.mass)
            

        self.save_snapshot()        
            
        # for plotting data
        e_in_array = np.array(e_in_array)
        g_in_array = np.array(g_in_array)
        o_in_array = np.array(o_in_array)
        e_out_array = np.array(e_out_array)
        g_out_array = np.array(g_out_array)
        o_out_array = np.array(o_out_array)
        i_relative_array = np.array(i_relative_array)

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
        self.plot_data.i_relative_array = i_relative_array
        self.plot_data.m1_array = m1_array
        self.plot_data.m2_array = m2_array
        self.plot_data.m3_array = m3_array
    #-------




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
    i_relative_array = triple.plot_data.i_relative_array
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
    plot_i_relative = figure.add_subplot(N_subplots,1,2)
    plot_e_in_g_in = figure.add_subplot(N_subplots,1,3)
    plot_a_in = figure.add_subplot(N_subplots,1,4)


    plot_e_in.plot(times_array_Myr,e_in_array)
    plot_e_in.set_xlim(0,t_max_Myr)
    plot_e_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_e_in.set_ylabel('$e_\mathrm{in}$')

    plot_i_relative.plot(times_array_Myr,i_relative_array*180.0/np.pi)
    plot_i_relative.set_xlim(0,t_max_Myr)
    plot_i_relative.set_ylim(0.9*min(i_relative_array*180.0/np.pi),1.1*max(i_relative_array*180.0/np.pi))
    plot_i_relative.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_relative.set_ylabel('$i_\mathrm{relative} ({}^\circ)$')

    plot_e_in_g_in.plot(np.cos(g_in_array),e_in_array)
    plot_e_in_g_in.set_xlabel('$\cos(g_\mathrm{in})$')
    plot_e_in_g_in.set_ylabel('$e_\mathrm{in}$')

    plot_a_in.plot(times_array_Myr,a_in_array_AU)
    plot_a_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_in.set_ylabel('$a_\mathrm{in}$')
    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)

    generic_plot_name = '_M'+str(m1_array[0]) + '_m'+str(m2_array[0]) +'_n'+str(m3_array[0]) + '_a'+str(a_in_array_AU[0]) + '_A'+str(a_out_array_AU[0]) + '_e'+str(e_in_array[0]) + '_E'+str(e_out_array[0]) + '_i'+str(i_relative_array[0]/np.pi*180.0) + '_g'+str(g_in_array[0]) + '_G'+str(g_out_array[0]) + '_o'+str(o_in_array[0]) + '_O'+str(o_out_array[0]) + '_t'+str(t_max_Myr) + '_maxdr'+str(maximum_radius_change_factor)
    plt.savefig('plots/orbit/TPS_inner_orbit'+generic_plot_name+'.pdf')
    plt.show()




    plt.plot(times_array_Myr,e_in_array)
    plt.plot(times_array_Myr,e_in_array, '.')
    plt.xlim(0,t_max_Myr)
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$e_\mathrm{in}$')
    plt.savefig('plots/orbit/e_in_time_'+generic_plot_name+'.pdf')
    plt.show()


#   wind a = ai * Mti/Mt
#    Mtot = m1_array+m2_array    
#    plt.plot(times_array_Myr,a_in_array_AU)
#    plt.plot(times_array_Myr,a_in_array_AU, '.')
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot)
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot, '.')
#    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*Mtot[:-1]/Mtot[1:])
#    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*Mtot[:-1]/Mtot[1:], '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$a_\mathrm{in}$')
#    plt.savefig('plots/orbit/semi_inner_wind'+generic_plot_name+'.pdf')
#    plt.show()
#    
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot/a_in_array_AU)
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot/a_in_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{in}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig('plots/orbit/semi_inner_rel_wind'+generic_plot_name+'.pdf')
#    plt.show()


#   cons mt a = ai * (m1i*m2i*/m1/m2)**2
#    plt.plot(times_array_Myr,a_in_array_AU, 'b-')
#    plt.plot(times_array_Myr,a_in_array_AU, 'b.')
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2, 'g-')
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2, 'g.')
#    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*(m1_array[:-1]*m2_array[:-1]/m1_array[1:]/m2_array[1:])**2, 'r-')
#    plt.plot(times_array_Myr[1:], a_in_array_AU[:-1]*(m1_array[:-1]*m2_array[:-1]/m1_array[1:]/m2_array[1:])**2, 'r.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$a_\mathrm{in}$')
#    plt.savefig('plots/orbit/semi_inner_cons_mt'+generic_plot_name+'.pdf')
#    plt.show()
#    
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2/a_in_array_AU)
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2/a_in_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{in}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig('plots/orbit/semi_inner_rel_cons_mt'+generic_plot_name+'.pdf')
#    plt.show()



#    dm = (m1_array[1:] - m1_array[:-1] )
#    dt = (times_array_Myr[1:] - times_array_Myr[:-1])
#    dmdt = (m1_array[1:] - m1_array[:-1] )/(times_array_Myr[1:] - times_array_Myr[:-1])

#    plt.plot(dmdt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:])
#    plt.plot(dmdt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:], '.')
#    plt.show()
#
#    plt.plot(dm, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:])
#    plt.plot(dm, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:], '.')
#    plt.show()
#
#    plt.plot(dt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:])
#    plt.plot(dt, a_in_array_AU[0]*Mtot[0]/Mtot[1:]/a_in_array_AU[1:], '.')
#    plt.show()
#
    
    
    
    
    
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
    plot_i_relative2 = figure.add_subplot(N_subplots,1,2)
    plot_e_out_g_out = figure.add_subplot(N_subplots,1,3)
    plot_a_out = figure.add_subplot(N_subplots,1,4)

#    times_array_Myr = triple.plot_data.times_array.value_in(units.Myr)
#    t_max_Myr = max(times_array_Myr)
#    i_relative_array = triple.plot_data.i_relative_array

    plot_e_out.plot(times_array_Myr,e_out_array)
    plot_e_out.set_xlim(0,t_max_Myr)
    plot_e_out.set_xlabel('$t/\mathrm{Myr}$')
    plot_e_out.set_ylabel('$e_\mathrm{out}$')

    plot_i_relative2.plot(times_array_Myr,i_relative_array*180.0/np.pi)
    plot_i_relative2.set_xlim(0,t_max_Myr)
    plot_i_relative2.set_ylim(0.9*min(i_relative_array*180.0/np.pi),1.1*max(i_relative_array*180.0/np.pi))
    plot_i_relative2.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_relative2.set_ylabel('$i_\mathrm{relative} ({}^\circ)$')

    plot_e_out_g_out.plot(np.cos(g_out_array),e_out_array)
    plot_e_out_g_out.set_xlabel('$\cos(g_\mathrm{out})$')
    plot_e_out_g_out.set_ylabel('$e_\mathrm{out}$')

    plot_a_out.plot(times_array_Myr,a_out_array_AU)
    plot_a_out.set_xlabel('$t/\mathrm{Myr}$')
    plot_a_out.set_ylabel('$a_\mathrm{out}$')

    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)
    plt.savefig('plots/orbit/TPS_outer_orbit'+generic_plot_name+'.pdf')
    plt.show()



#    plt.plot(times_array_Myr,e_out_array)
#    plt.xlim(0,t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$e_\mathrm{out}$')
#    plt.savefig('plots/orbit/e_out_time_'+generic_plot_name+'.pdf')
#    plt.show()


#    Mtott = m1_array+m2_array+m3_array    
#    plt.plot(times_array_Myr,a_out_array_AU)
#    plt.plot(times_array_Myr,a_out_array_AU, '.')
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott)
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott, '.')
#    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*Mtott[:-1]/Mtott[1:])
#    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*Mtott[:-1]/Mtott[1:], '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$a_\mathrm{out}$')
#    plt.savefig('plots/orbit/semi_outer_wind'+generic_plot_name+'.pdf')
#    plt.show()
#
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott/a_out_array_AU)
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott/a_out_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{out}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig('plots/orbit/semi_outer_rel_wind'+generic_plot_name+'.pdf')
#    plt.show()
#
#    m_in_array = m1_array+m2_array
#    plt.plot(times_array_Myr,a_out_array_AU, 'b-')
#    plt.plot(times_array_Myr,a_out_array_AU, 'b.')
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2, 'g-')
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2, 'g.')
#    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*(m_in_array[:-1]*m3_array[:-1]/m_in_array[1:]/m3_array[1:])**2, 'r-')
#    plt.plot(times_array_Myr[1:], a_out_array_AU[:-1]*(m_in_array[:-1]*m3_array[:-1]/m_in_array[1:]/m3_array[1:])**2, 'r.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$a_\mathrm{out}$')
#    plt.savefig('plots/orbit/semi_outer_cons_mt'+generic_plot_name+'.pdf')
#    plt.show()
#
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2/a_out_array_AU)
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2/a_out_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{out}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig('plots/orbit/semi_outer_rel_cons_mt'+generic_plot_name+'.pdf')
#    plt.show()
#

#    dm = (m3_array[1:] - m3_array[:-1] )
#    dmdt = (m3_array[1:] - m3_array[:-1] )/(times_array_Myr[1:] - times_array_Myr[:-1])

#    plt.plot(dmdt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:])
#    plt.plot(dmdt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:], '.')
#    plt.show()
#
#    plt.plot(dm, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:])
#    plt.plot(dm, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:], '.')
#    plt.show()
#
#    plt.plot(dt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:])
#    plt.plot(dt, a_out_array_AU[0]*Mtott[0]/Mtott[1:]/a_out_array_AU[1:], '.')
#    plt.show()



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
            relative_inclination= 80.0*np.pi/180.0,
            inner_argument_of_pericenter= 0.1, outer_argument_of_pericenter= 0.5,
            inner_longitude_of_ascending_node= 0.0, outer_longitude_of_ascending_node= 0.0,
            metallicity= 0.02,
            tend= 5.0 |units.Myr, number = 0):

    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    inner_eccentricity = float(inner_eccentricity)
    outer_eccentricity = float(outer_eccentricity)
    relative_inclination = float(relative_inclination)
    inner_argument_of_pericenter = float(inner_argument_of_pericenter)
    outer_argument_of_pericenter = float(outer_argument_of_pericenter)
    inner_longitude_of_ascending_node = float(inner_longitude_of_ascending_node)
    outer_longitude_of_ascending_node = float(outer_longitude_of_ascending_node)

    triple_class_object = Triple_Class(inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            relative_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node,
            metallicity, tend, number)

    triple_class_object.evolve_model()
    plot_function(triple_class_object)
#    triple_class_object.print_stellar_system()
    return triple_class_object
#-----

#-----
#for running triple.py from the commandline
def parse_arguments():
    from amuse.units.optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-M", unit=units.MSun, 
                      dest="inner_primary_mass", type="float", default = 1.3|units.MSun,
                      help="inner primary mass [%default]")
    parser.add_option("-m",  unit=units.MSun, 
                      dest="inner_secondary_mass", type="float", default = 0.5|units.MSun,
                      help="inner secondary mass [%default]")
    parser.add_option("-l",  unit=units.MSun, 
                      dest="outer_mass", type="float", default = 0.5|units.MSun,
                      help="outer mass [%default]")

    parser.add_option("-A", unit=units.RSun,
                      dest="inner_semimajor_axis", type="float", 
                      default = 200.0 |units.RSun,
                      help="inner semi major axis [%default]")
    parser.add_option("-a", unit=units.RSun,
                      dest="outer_semimajor_axis", type="float", 
                      default = 20000.0 |units.RSun,
                      help="outer semi major axis [%default]")
    parser.add_option("-E",
                      dest="inner_eccentricity", type="float", default = 0.1,
                      help="inner eccentricity [%default]")
    parser.add_option("-e",
                      dest="outer_eccentricity", type="float", default = 0.5,
                      help="outer eccentricity [%default]")
    parser.add_option("-i","-I",
                      dest="relative_inclination", type="float", default = 80.0*np.pi/180.0,
                      help="relative inclination [rad] [%default]")
    parser.add_option("-G",
                      dest="inner_argument_of_pericenter", type="float", default = 0.1,
                      help="inner argument of pericenter [rad] [%default]")
    parser.add_option("-g",
                      dest="outer_argument_of_pericenter", type="float", default = 0.5,
                      help="outer argument of pericenter [rad] [%default]")
    parser.add_option("-O",
                      dest="inner_longitude_of_ascending_node", type="float", default = 0.0,
                      help="inner longitude of ascending node [rad] [%default]")
    parser.add_option("-o",
                      dest="outer_longitude_of_ascending_node", type="float", default = 0.0,
                      help="outer longitude of ascending node [rad] [%default]")
    parser.add_option("-z", dest="metallicity", type="float", default = 0.02,
                      help="metallicity [%default] %unit")
    parser.add_option("-t", "-T", unit=units.Myr, 
                      dest="tend", type="float", default = 5.0 |units.Myr,
                      help="end time [%default] %unit")
    parser.add_option("-N", dest="number", type="int", default = 0,
                      help="number of system [%default]")

    options, args = parser.parse_args()
    return options.__dict__





if __name__ == '__main__':
    options = parse_arguments()

    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")


    triple_class_object = Triple_Class(**options)
    triple_class_object.evolve_model()
#    triple_class_object.print_stellar_system()
#    plot_function(triple_class_object)


