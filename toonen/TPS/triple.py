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

import os, sys
sys.path.insert(1, os.environ.get('AMUSE_DIR')+'/sandbox/hamers/TPS/code')
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
REPORT_DT = False 

#stopping conditions
#for the moment these are global variables, and not interface parameters, as we cannot simulate them
stop_at_merger = True # as implementation is missing, e.g. secular code not adjusted
stop_at_disintegrated = True # as implementation is missing, e.g. secular code not adjusted
stop_at_triple_mass_transfer = True # as implementation is missing 
stop_at_collision = True # as implementation is missing
stop_at_dynamical_instability = True # this should always be true!

stop_at_mass_transfer = False
no_stellar_evolution = False

file_name = "triple.hdf" #"triple.txt"
file_type = "hdf5" #"txt"

#constants
time_step_factor_stable_mt = 0.01 #1% mass loss during mass transfer
# lowering this to 0.005 makes the code twice as slow
# 0.01 -> error in the semi-major axis of about 0.5%
maximum_wind_mass_loss_factor = 0.01 
error_dm = 0.05
#maximum_radius_change_factor = 0.005
error_dr = 0.01
minimum_time_step = 1.e-9 |units.Myr
min_mass = 0.08 |units.MSun # for stars
max_mass = 100 |units.MSun
maximum_time_step_factor = 100.
#Rl_fraction = 0.9#1.0-10.*error_dr # ratio or star radius over Roche lobe at which time step is decreased
                              # radius grows maximally by error_dr
time_step_factor_kozai = 0.025 # 0.2*0.1, 0.2-> for error in kozai timescale, 0.1-> 10 steps per cycle



stellar_types_SN_remnants = [13,14]|units.stellar_type # remnant types created through a supernova
stellar_types_remnants = [7,8,9,10,11,12,13,14]|units.stellar_type


class Triple_Class:
    #-------
    #setup stellar system
    def __init__(self, inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            relative_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node,
            metallicity, tend, number, maximum_radius_change_factor, tidal_terms):      
            
        inner_eccentricity, outer_eccentricity = self.test_initial_parameters(inner_primary_mass, inner_secondary_mass, outer_mass,
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
        
        self.first_contact = False 
        self.instantaneous_evolution = False # no secular evolution        
        self.tend = tend #...
        self.time = 0.0|units.yr
        self.previous_time = 0.0|units.yr
        self.maximum_radius_change_factor = maximum_radius_change_factor

        self.triple = bins[1]
        self.triple.relative_inclination = relative_inclination 
        self.triple.is_star = False
        self.triple.dynamically_stable = True 
        self.triple.number = number 
        
        self.setup_stellar_code(metallicity, stars)
        self.setup_secular_code(self.triple.as_set(), tidal_terms)

        self.triple.kozai_type = self.get_kozai_type()
        self.update_previous_stellar_parameters()
        self.update_stellar_parameters() 
        self.update_time_derivative_of_radius()


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
#        self.stellar_code = SeBa(redirection='none')

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

      
    def setup_secular_code(self, triple_set, tidal_terms):
        self.secular_code = SecularTriple()
#        self.secular_code = SecularTriple(redirection='none')
        self.secular_code.triples.add_particles(triple_set)
        self.secular_code.parameters.equations_of_motion_specification = 0
        self.secular_code.parameters.include_quadrupole_terms = True
        self.secular_code.parameters.include_octupole_terms = True        
        self.secular_code.parameters.include_inner_wind_terms = True
        self.secular_code.parameters.include_outer_wind_terms = True
        self.secular_code.parameters.include_inner_RLOF_terms = True
        self.secular_code.parameters.include_outer_RLOF_terms = True
        self.secular_code.parameters.include_spin_radius_mass_coupling_terms = True
        self.secular_code.parameters.include_magnetic_braking_terms = False # not tested

        self.secular_code.parameters.include_inner_tidal_terms = tidal_terms
        self.secular_code.parameters.include_outer_tidal_terms = tidal_terms
        
        self.secular_code.parameters.include_1PN_inner_terms = False
        self.secular_code.parameters.include_1PN_outer_terms = False
        self.secular_code.parameters.include_1PN_inner_outer_terms = False ### warning: probably broken
        self.secular_code.parameters.include_25PN_inner_terms = False
        self.secular_code.parameters.include_25PN_outer_terms = False

        self.secular_code.parameters.check_for_dynamical_stability = True
        self.secular_code.parameters.check_for_inner_collision = True
        self.secular_code.parameters.check_for_outer_collision = True

#   gives problems when RLOF begins or starts during timestep in secular code        
#        self.secular_code.parameters.check_for_inner_RLOF = True 
#        self.secular_code.parameters.check_for_outer_RLOF = True 

         # accuracy of secular code
#        self.secular_code.parameters.input_precision = 1.0e-10#1.0e-5
#        self.secular_code.parameters.relative_tolerance = 1.0e-10
        self.secular_code.parameters.include_linear_mass_change = True #needed for Jspin conservation
        self.secular_code.parameters.include_linear_radius_change = True #needed for Jspin conservation

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
        if (inner_eccentricity < minimum_eccentricity):
            inner_eccentricity = minimum_eccentricity
        if (outer_eccentricity < minimum_eccentricity):
            outer_eccentricity = minimum_eccentricity
    
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
            
        return inner_eccentricity, outer_eccentricity            
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
            if self.time == quantities.zero: #initialization
               stellar_system.previous_time_derivative_of_radius = 0.0 | units.RSun/units.yr
            else:
               stellar_system.previous_time_derivative_of_radius = stellar_system.time_derivative_of_radius
        else:
            self.update_previous_stellar_parameters(stellar_system.child1)        
            self.update_previous_stellar_parameters(stellar_system.child2)
            stellar_system.previous_mass = self.get_mass(stellar_system) 
            
            if self.is_triple(stellar_system):
                stellar_system.previous_kozai_type = stellar_system.kozai_type
    #-------

    #-------
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
            stellar_system.gyration_radius = stellar_system.gyration_radius_sq**0.5     
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
            if self.is_triple(stellar_system):
                stellar_system.kozai_type = self.get_kozai_type()
            
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
    
    

    def has_kozai_type_changed(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if self.is_triple(stellar_system):
            if stellar_system.kozai_type != stellar_system.previous_kozai_type:
                return True
        else: #binaries and single stars do not have a kozai timescale
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
        if not bs.is_star:
            Porb = 2*np.pi * np.sqrt(bs.semimajor_axis**3/constants.G / self.get_mass(bs))
            return Porb
        else:
            print 'orbital_period: single star does not have a period'
            exit(-1)

    def orbital_angular_momentum(self, bs):
        if not bs.is_star:
            M = self.get_mass(bs.child1)
            m = self.get_mass(bs.child2)
            a = bs.semimajor_axis
            e = bs.eccentricity
            J = M*m * np.sqrt(constants.G*a*(1-e**2)/(M+m))
        
            if REPORT_BINARY_EVOLUTION:
                print 'Jorb:', M, m, a, e, J
        
            return J
        else:
            print 'orbital_angular_momentum: single star does not have an orbit'
            exit(-1)
    
    def spin_angular_momentum(self, ss):
        if ss.is_star:
            moment_of_inertia = ss.gyration_radius**2 * ss.mass * ss.radius**2
            Jstar = moment_of_inertia * ss.spin_angular_frequency
            return Jstar            
        else:
            print 'spin_angular_momentum: structure stellar system unknown'        
            exit(2)
            
    def apsidal_motion_constant(self, star):
        if star.stellar_type in [13, 14]|units.stellar_type: #ns
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
            
    def get_kozai_type(self):
        if self.is_triple():
            t_kozai = self.kozai_timescale()
            t_nuc = self.get_max_nuclear_evolution_timescale_of_system()
            if t_kozai < t_nuc:
                return True
            else:
                return False
        else:
           print 'Kozai type needs triple system'
           exit(1)   

    def get_max_nuclear_evolution_timescale_of_system(self, stellar_system = None):
        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            return nuclear_evolution_timescale(stellar_system)
        else:
            t1 = self.get_max_nuclear_evolution_timescale_of_system(stellar_system.child1)        
            t2 = self.get_max_nuclear_evolution_timescale_of_system(stellar_system.child2)
            return max(t1, t2)


    def check_for_RLOF(self):
        if self.triple.is_star:
            return
        elif self.is_binary():
            Rl1 = roche_radius(self, self.child1)
            Rl2 = roche_radius(self, self.child2)
            if REPORT_TRIPLE_EVOLUTION:
                print 'Roche lobe radii:', Rl1, Rl2
                print 'Stellar radii:', self.triple.child1.radius, self.triple.child2.radius
            
            first_RLOF = []
            if self.triple.child1.radius >= Rl1:
                if not self.triple.child1.is_donor:
                    first_RLOF.append(True)
                self.triple.child1.is_donor = True
            else:                                                         
               self.triple.child1.is_donor = False

            if self.triple.child2.radius >= Rl2:
                if not self.triple.child2.is_donor:
                    first_RLOF.append(True)
                self.triple.child2.is_donor = True
            else:                                                         
               self.triple.child2.is_donor = False
               
            self.first_contact = False
            if any(first_RLOF):
                self.first_contact = True
               
 
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
                print 'binary Roche lobe radii:', roche_radius(bin, bin.child1, self), roche_radius(bin, bin.child2, self)

            first_RLOF = []
            if bin.child1.radius >= Rl1:
                if not bin.child1.is_donor:
                    first_RLOF.append(True)
                bin.child1.is_donor = True
            else:
                bin.child1.is_donor = False
            
            if bin.child2.radius >= Rl2:
                if not bin.child2.is_donor:
                    first_RLOF.append(True)
                bin.child2.is_donor = True
            else:
                bin.child2.is_donor = False
            
            if star.radius >= Rl3:
                if not star.is_donor:
                    first_RLOF.append(True)
                star.is_donor = True
            else:
                star.is_donor = False

            self.first_contact = False
            if any(first_RLOF):
                self.first_contact = True
                
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

            if REPORT_TRIPLE_EVOLUTION:
                print 'dt radius change:', stellar_system.time_derivative_of_radius, stellar_system.previous_time_derivative_of_radius 
                print stellar_system.stellar_type, stellar_system.previous_stellar_type
                print stellar_system.radius, stellar_system.previous_radius
                print stellar_system.mass, stellar_system.previous_mass            
            
            if stellar_system.wind_mass_loss_rate * -1. > quantities.zero:
                dt = maximum_wind_mass_loss_factor*stellar_system.mass / stellar_system.wind_mass_loss_rate*-1.
            if REPORT_DT:
                print "Dt_wind_star = ", dt
            return dt 
        else:
            dt1 = self.determine_time_step_wind(stellar_system.child1)        
            dt2 = self.determine_time_step_wind(stellar_system.child2)
            if REPORT_DT:
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
            if stellar_system.time_derivative_of_radius != quantities.zero:
                growth_factor = 1. # Rdot < Rdot_prev
                if stellar_system.time_derivative_of_radius > stellar_system.previous_time_derivative_of_radius:
                    if stellar_system.previous_time_derivative_of_radius == quantities.zero:
                        growth_factor = 0.1
                    elif stellar_system.previous_time_derivative_of_radius * stellar_system.time_derivative_of_radius < quantities.zero:
                        growth_factor = 0.01
                    elif stellar_system.time_derivative_of_radius > quantities.zero:
                        growth_factor = stellar_system.previous_time_derivative_of_radius/stellar_system.time_derivative_of_radius 
                    else:
                        growth_factor = stellar_system.time_derivative_of_radius/stellar_system.previous_time_derivative_of_radius 


#                growth_factor = 0.1
#                if stellar_system.previous_time_derivative_of_radius != quantities.zero:
#                    growth_factor = stellar_system.previous_time_derivative_of_radius/stellar_system.time_derivative_of_radius
#                    if abs(growth_factor) > 1: 
#                        growth_factor = 1./growth_factor 
#                    if stellar_system.previous_time_derivative_of_radius * stellar_system.time_derivative_of_radius < quantities.zero:
#                        growth_factor = growth_factor * 0.01                        
#                                   
#                print 'dt radius change:', stellar_system.time_derivative_of_radius, stellar_system.previous_time_derivative_of_radius 
#                print stellar_system.stellar_type, stellar_system.previous_stellar_type
#                print growth_factor, stellar_system.radius
#                print stellar_system.mass, stellar_system.previous_mass

                dt = abs(growth_factor * self.maximum_radius_change_factor*stellar_system.radius / stellar_system.time_derivative_of_radius)


            if REPORT_DT:
                print "Dt_radius_change_star = ", dt
            return dt 
        else:
            dt1 = self.determine_time_step_radius_change(stellar_system.child1)        
            dt2 = self.determine_time_step_radius_change(stellar_system.child2)
            if REPORT_DT:
                print "Dt_radius_change_binary = ", dt1, dt2
            return min(dt1, dt2) 

    def determine_time_step_kozai(self):
    #note: returned value can be inf when the system is a binary or single star            
        if self.is_triple():
            dt = self.kozai_timescale()*time_step_factor_kozai
            if REPORT_DT:
                print "Dt_kozai = ", dt
        else:
            dt = np.inf |units.Myr
        
        P_out = self.orbital_period(self.triple)
        P_in = self.orbital_period(self.triple.child2)        
        return dt

 
    def determine_time_step_stable_mt(self, stellar_system = None):
    #note: returned value can be inf when the mass_transfer_rate = 0 
        if stellar_system == None:
            stellar_system = self.triple
            
        if stellar_system.is_star:
            dt = np.inf |units.Myr
            if stellar_system.is_donor:
                print 'dt_mt',time_step_factor_stable_mt, stellar_system.mass, stellar_system.parent.mass_transfer_rate
                dt = abs(time_step_factor_stable_mt*stellar_system.mass/stellar_system.parent.mass_transfer_rate)
            if REPORT_DT:
                print "Dt_mt_star = ", dt
            return dt
        else:
            if self.is_binary(stellar_system) and stellar_system.child1.is_donor and stellar_system.child2.is_donor:
                #should have been taken care of in determine_time_step()
                print 'determine_time_step_stable_mt: contact system'        
                exit(1)

            dt1 = self.determine_time_step_stable_mt(stellar_system.child1)        
            dt2 = self.determine_time_step_stable_mt(stellar_system.child2)
            if REPORT_DT:
                print "Dt_mt_binary = ", dt1, dt2
            return min(dt1, dt2) 
         
#    def close_to_RLOF(self):
#
#        if self.triple.is_star:
#            return False, 0        
#        elif self.is_binary():
#            Rl1 = roche_radius(self, self.child1)
#            Rl2 = roche_radius(self, self.child2)
#            if self.triple.child1.radius >= Rl1 or self.triple.child2.radius >= Rl2:
#                print 'close to rlof already rlof'
#                exit(0)
#            elif self.triple.child1.radius >= Rl_fraction*Rl1 and self.triple.child2.radius >= Rl_fraction*Rl2:
#                ratio_rad_rlof = maximum(self.triple.child1.radius/Rl1, self.triple.child2.radius/Rl2)
#                print 'dt_close_to_mt', time_step_factor_stable_mt, self.triple.child1.mass, self.triple.child2.mass, self.triple.mass_transfer_rate
#                dt = abs(time_step_factor_stable_mt*min(self.triple.child1.mass, self.triple.child2.mass)/self.triple.mass_transfer_rate)
#                return True, ratio_rad_rlof, dt  
#            elif self.triple.child1.radius >= Rl_fraction*Rl1:
#                ratio_rad_rlof = self.triple.child1.radius/Rl1
#                print 'dt_close_to_mt', time_step_factor_stable_mt, self.triple.child1.mass, self.triple.mass_transfer_rate
#                dt = abs(time_step_factor_stable_mt*self.triple.child1.mass/self.triple.mass_transfer_rate)
#                return True, ratio_rad_rlof, dt  
#            elif self.triple.child2.radius >= Rl_fraction*Rl2:
#                ratio_rad_rlof = self.triple.child2.radius/Rl2
#                print 'dt_close_to_mt', time_step_factor_stable_mt, self.triple.child2.mass, self.triple.mass_transfer_rate
#                dt = abs(time_step_factor_stable_mt*self.triple.child2.mass/self.triple.mass_transfer_rate)
#                return True, ratio_rad_rlof, dt              
#            else:
#                return False, 0, 0    
#        elif self.is_triple():
#            if self.triple.child1.is_star:
#                star = self.triple.child1
#                bin = self.triple.child2
#            else:
#                star = self.triple.child2
#                bin = self.triple.child1
#
#            #assumping secular code always returns inner binary first
#            Rl1, Rl2, Rl3 = self.secular_code.give_roche_radii(self.triple)
#            if bin.child1.radius >= Rl1 or bin.child2.radius >= Rl2 or star.radius >= Rl3:
#                print 'close to rlof already rlof'
#                exit(0)
#
#            dt=[]
#            ratio_rad_rlof=[]
#            if star.radius >= Rl_fraction*Rl3:
#                ratio_rad_rlof.append(star.radius / Rl3)
#                print 'dt_close_to_mt', time_step_factor_stable_mt, star.mass, star.parent.mass_transfer_rate
#                dt.append(abs(time_step_factor_stable_mt*star.mass/star.parent.mass_transfer_rate))
#            if bin.child1.radius >= Rl_fraction*Rl1:
#                ratio_rad_rlof.append(bin.child1.radius / Rl1)
#                print 'dt_close_to_mt', time_step_factor_stable_mt, bin.child1.mass, bin.mass_transfer_rate
#                dt.append(abs(time_step_factor_stable_mt*bin.child1.mass/bin.mass_transfer_rate))
#            if bin.child2.radius >= Rl_fraction*Rl2:
#                ratio_rad_rlof.append(bin.child2.radius / Rl2)
#                print 'dt_close_to_mt', time_step_factor_stable_mt, bin.child2.mass, bin.mass_transfer_rate
#                dt.append(abs(time_step_factor_stable_mt*bin.child2.mass/bin.mass_transfer_rate))
#                    
#            if len(dt) > 0:
#                print 'close to rlof', ratio_rad_rlof, dt
#                return True, max(ratio_rad_rlof), min(dt)
#            else:
#                return False, 0, 0
#        else:
#            print 'close_to_RLOF: structure stellar system unknown'        
#            exit(2)    
    
        
    def determine_time_step(self):         
        if REPORT_DT:
            print "Dt = ", self.stellar_code.particles.time_step, self.tend/100.0

#        during unstable mass transfer, contact_system or other instantaneous interactions: no step in time
        if self.has_donor() and (not self.is_system_stable() or self.has_contact_system()):
            return minimum_time_step#0.0|units.yr -> problem for time_derivative_of_radius         
                       
        #maximum time_step            
        time_step_max = self.tend - self.time
        
        # time_step of stellar evolution
        time_step_stellar_code = self.stellar_code.particles.time_step
                    
        # small time_step during heavy wind mass losses
        time_step_wind =  self.determine_time_step_wind()

        # small time_step during phases of fast growth (note: not yet during fast shrinkage)
        time_step_radius_change = self.determine_time_step_radius_change()             
        
        #small time_step when kozai cycles affect the orbit significantly
        time_step_kozai = self.determine_time_step_kozai()
        
        time_step = min(time_step_kozai, min(time_step_radius_change, min(time_step_wind, min( min(time_step_stellar_code), time_step_max))))    
        if REPORT_DT or REPORT_TRIPLE_EVOLUTION:
            print 'time:', time_step_max, time_step_stellar_code, time_step_wind, time_step_radius_change, time_step_kozai, time_step


                
        #during stable mass transfer   
        if self.has_donor():
            print time_step, self.determine_time_step_stable_mt()
            time_step = min(time_step, self.determine_time_step_stable_mt())
#        else:
#            #during run-up towards mass transfer 
#            close_to_RLOF_bool, ratio_rad_rlof, time_step_close_to_mt = self.close_to_RLOF()     
#            if close_to_RLOF_bool:
#                if ratio_rad_rlof >= 1.0:
#                    print 'ratio_rad_rlof >1?'
#                    exit(0)
#                else:
#                    print 'close to rlof', time_step_close_to_mt, time_step, ratio_rad_rlof
#                    time_step = min(time_step, time_step_close_to_mt)

           
        if self.time == quantities.zero:
            print 'First timestep - outer period'
            #initialization (e.g. time_derivative_of_radius)
            P_out = self.orbital_period(self.triple) #period outer binary 
            # do not take 0.1*P_in -> resonance -> large error
            time_step = min(min(P_out, time_step), 1.|units.yr)
        else:
            previous_time_step = self.time - self.previous_time
#            print previous_time_step, maximum_time_step_factor*previous_time_step, time_step
            time_step = min(time_step, maximum_time_step_factor*previous_time_step)  
                    

#        if time_step < minimum_time_step:
#            print 'error small time_step'
#            print time_step_max, time_step_stellar_code, time_step_wind, time_step_radius_change, time_step
        time_step = max(time_step, minimum_time_step)  
        
        return time_step
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

    
    def safety_check_time_step(self, stellar_system = None):

        if stellar_system == None:
            stellar_system = self.triple

        if stellar_system.is_star:
            dm = (stellar_system.previous_mass - stellar_system.mass)/stellar_system.mass
            if REPORT_TRIPLE_EVOLUTION:
                print 'wind mass loss rate:', stellar_system.wind_mass_loss_rate,
                print 'relative wind mass losses:', dm

            if (dm > error_dm) and not (stellar_system.stellar_type != stellar_system.previous_stellar_type and stellar_system.stellar_type in stellar_types_SN_remnants):
                print 'Change in mass in a single time_step larger then', error_dm
                print dm, stellar_system.mass, stellar_system.stellar_type
                exit(1)

            if self.secular_code.parameters.include_inner_tidal_terms or self.secular_code.parameters.include_outer_tidal_terms:
                dr = (stellar_system.radius - stellar_system.previous_radius)/stellar_system.radius
    
                if REPORT_TRIPLE_EVOLUTION:    
                    print 'change in radius over time:',  stellar_system.time_derivative_of_radius,
                    print 'relative change in radius:', dr
    
                if (dr > error_dr) and not (stellar_system.stellar_type != stellar_system.previous_stellar_type and stellar_system.stellar_type in stellar_types_remnants):
                    print 'Change in radius in a single time_step larger then', error_dr
                    print dr, stellar_system.time_derivative_of_radius, stellar_system.mass, stellar_system.previous_mass, stellar_system.stellar_type, stellar_system.previous_stellar_type, self.has_stellar_type_changed(stellar_system)
                    exit(1)

        else:
            self.safety_check_time_step(stellar_system.child1)        
            self.safety_check_time_step(stellar_system.child2)
            
    
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
        spin1_array = []
        spin2_array = []
        spin3_array = []
        r1_array = []
        r2_array = []
        r3_array = []
        gy1_array = []
        gy2_array = []
        gy3_array = []
    
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
        spin1_array.append(self.triple.child2.child1.spin_angular_frequency.value_in(1./units.Myr))
        spin2_array.append(self.triple.child2.child2.spin_angular_frequency.value_in(1./units.Myr))
        spin3_array.append(self.triple.child1.spin_angular_frequency.value_in(1./units.Myr))
        r1_array.append(self.triple.child2.child1.radius.value_in(units.RSun))
        r2_array.append(self.triple.child2.child2.radius.value_in(units.RSun))
        r3_array.append(self.triple.child1.radius.value_in(units.RSun))
        gy1_array.append(self.triple.child2.child1.gyration_radius.number)
        gy2_array.append(self.triple.child2.child2.gyration_radius.number)
        gy3_array.append(self.triple.child1.gyration_radius.number)
        
        if REPORT_TRIPLE_EVOLUTION:
            print 'kozai timescale:', self.kozai_timescale(), self.triple.kozai_type, self.tend    
        print 'kozai timescale:', self.kozai_timescale(), self.triple.kozai_type, self.tend    
        self.determine_mass_transfer_timescale()
        self.save_snapshot()    
        while self.time<self.tend:
            if self.has_stellar_type_changed() or self.has_kozai_type_changed():
                self.save_snapshot()        
   
            dt = self.determine_time_step()  
            if not no_stellar_evolution: 
                self.update_previous_stellar_parameters()
    
            self.time += dt            
            
            if REPORT_TRIPLE_EVOLUTION:
                print '\n\ntime:', self.time, self.has_donor()            

            #do stellar evolution 
            if not no_stellar_evolution: 
                if REPORT_TRIPLE_EVOLUTION:
                    print 'Stellar evolution'

                self.stellar_code.evolve_model(self.time)
                self.channel_from_stellar.copy()
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
            
            if stop_at_mass_transfer and self.has_donor():
                print "Mass transfer at time/Myr = ",self.time.value_in(units.Myr)
                self.save_snapshot()        
                break             
            if stop_at_triple_mass_transfer and self.has_triple_mass_transfer():
                print 'Mass transfer in outer binary of triple at time/Myr = ",self.time.value_in(units.Myr)'
                #it's possible that there is mass transfer in the inner and outer binary
                print self.triple.child2.bin_type
                print self.triple.bin_type
                print self.has_triple_mass_transfer()
                break                                    
                                    
            #do stellar interaction
            if REPORT_TRIPLE_EVOLUTION:
                print 'Stellar interaction'
            self.determine_mass_transfer_timescale()
            self.resolve_stellar_interaction()
            self.update_time_derivative_of_radius() # includes radius change from wind and mass transfer

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
                print 'Merger at time/Myr = ",self.time.value_in(units.Myr)'
                break    
            if stop_at_disintegrated and self.has_disintegrated():
                print 'Disintegration of system at time/Myr = ",self.time.value_in(units.Myr)'
                break
            
            if self.instantaneous_evolution == False and self.first_contact == False: 
                # do secular evolution
                if REPORT_TRIPLE_EVOLUTION:
                    print 'Secular evolution'
                self.channel_to_secular.copy()
                self.safety_check_time_step()
                if not self.is_triple:# e.g. binaries
                    print 'Secular code disabled'
                    exit(1)

                # if mass transfer should only take place for a fraction of the timestep
                # e.g. at the end of mass transfer when the envelope is thin
                if self.triple.child2.part_dt_mt < 1: # inner binary, see function determine_partial_time_step_stable_mass_transfer
                    full_dt = self.time - self.previous_time
                    self.secular_code.evolve_model(self.previous_time + full_dt * self.triple.child2.part_dt_mt)
                    if self.triple.error_flag_secular < 0:
                        print "Error in secular code at time/Myr = ",self.time.value_in(units.Myr)
                        print self.triple.error_flag_secular
                        self.save_snapshot()        
                        break
                    self.channel_from_secular.copy()

                    self.check_for_RLOF() 
                    if self.has_donor():
                        self.check_for_RLOF()
                        print self.has_donor()
                        self.channel_from_secular.copy()
                        self.check_for_RLOF()
                        print 'After partial timestep the system should be detached...'
                        print self.has_donor()
                        print self.triple.child1.mass, self.triple.child2.child1.mass, self.triple.child2.child2.mass
                        break
                    # not necessary because secular code reset is_donor and therefore the mass transfer rate is not used if the system is detached    
#                    self.triple.child2.mass_transfer_rate =  0.0 | units.MSun/units.yr 
                        
                    self.secular_code.evolve_model(self.time)
                    self.triple.child2.part_dt_mt = 1.                    
                    
                else:
                    self.secular_code.evolve_model(self.time)
                    
                if REPORT_TRIPLE_EVOLUTION:
                    print 'Secular evolution finished'
                    

                if self.time - self.secular_code.model_time < -1*numerical_error|units.Myr:
                    print 'triple time < sec time: should not be possible', self.time, self.secular_code.model_time
                    exit(1)

                elif self.time - self.secular_code.model_time > numerical_error|units.Myr:
                 #When the secular code discovers RLOF or the end of RLOF, the orbital simulation stops. 
                 #The stellar evolution has evolved to far in time. For now this is not compensated. 
                    self.secular_code.model_time = self.time                     

                if stop_at_mass_transfer and self.has_donor():
                    print "Mass transfer at time/Myr = ",self.time.value_in(units.Myr)
                    self.save_snapshot()        
                    break             
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
                if self.triple.error_flag_secular < 0:
                    print "Error in secular code at time/Myr = ",self.time.value_in(units.Myr)
                    print self.triple.error_flag_secular
                    self.save_snapshot()        
                    break
                   
            else:
                print 'skip secular', self.instantaneous_evolution, self.first_contact        
                print self.instantaneous_evolution, self.first_contact
                self.secular_code.model_time = self.time
                self.instantaneous_evolution = False
                
            #should also do safety check time_step here -> only make sure that mass loss from stable mass transfer is not too large -> determine_time_step_mt
            self.check_for_RLOF()

            
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
            spin1_array.append(self.triple.child2.child1.spin_angular_frequency.value_in(1./units.Myr))
            spin2_array.append(self.triple.child2.child2.spin_angular_frequency.value_in(1./units.Myr))
            spin3_array.append(self.triple.child1.spin_angular_frequency.value_in(1./units.Myr))
            r1_array.append(self.triple.child2.child1.radius.value_in(units.RSun))
            r2_array.append(self.triple.child2.child2.radius.value_in(units.RSun))
            r3_array.append(self.triple.child1.radius.value_in(units.RSun))
            gy1_array.append(self.triple.child2.child1.gyration_radius.number)
            gy2_array.append(self.triple.child2.child2.gyration_radius.number)
            gy3_array.append(self.triple.child1.gyration_radius.number)
            

        self.save_snapshot()        
            
        # for plotting data
        e_in_array = np.array(e_in_array)
        g_in_array = np.array(g_in_array)
        o_in_array = np.array(o_in_array)
        e_out_array = np.array(e_out_array)
        g_out_array = np.array(g_out_array)
        o_out_array = np.array(o_out_array)
        i_relative_array = np.array(i_relative_array)
        spin1_array = np.array(spin1_array)
        spin2_array = np.array(spin2_array)
        spin3_array = np.array(spin3_array)
        r1_array = np.array(r1_array)
        r2_array = np.array(r2_array)
        r3_array = np.array(r3_array)
        gy1_array = np.array(gy1_array)
        gy2_array = np.array(gy2_array)
        gy3_array = np.array(gy3_array)

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
        self.plot_data.spin1_array = spin1_array
        self.plot_data.spin2_array = spin2_array
        self.plot_data.spin3_array = spin3_array
        self.plot_data.r1_array = r1_array
        self.plot_data.r2_array = r2_array
        self.plot_data.r3_array = r3_array
        self.plot_data.gy1_array = gy1_array
        self.plot_data.gy2_array = gy2_array
        self.plot_data.gy3_array = gy3_array
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
    spin1_array = triple.plot_data.spin1_array
    spin2_array = triple.plot_data.spin2_array
    spin3_array = triple.plot_data.spin3_array
    r1_array = triple.plot_data.r1_array
    r2_array = triple.plot_data.r2_array
    r3_array = triple.plot_data.r3_array
    gy1_array = triple.plot_data.gy1_array
    gy2_array = triple.plot_data.gy2_array
    gy3_array = triple.plot_data.gy3_array
    
    
    generic_name = 'M'+str(m1_array[0]) + '_m'+str(m2_array[0]) +'_n'+str(m3_array[0]) + '_a'+str(a_in_array_AU[0]) + '_A'+str(a_out_array_AU[0]) + '_e'+str(e_in_array[0]) + '_E'+str(e_out_array[0]) + '_i'+str(i_relative_array[0]/np.pi*180.0) + '_g'+str(g_in_array[0]) + '_G'+str(g_out_array[0]) + '_o'+str(o_in_array[0]) + '_O'+str(o_out_array[0]) + '_t'+str(t_max_Myr) + '_maxdr'+str(triple.maximum_radius_change_factor)+'_edr'+str(error_dr)
#    f = open(generic_name+'.txt','w')
#    for i in range(len(times_array_Myr)):
#        f.write( str(times_array_Myr[i]) + '\t'+str(e_in_array[i]) + '\t'+ str(a_in_array_AU[i]) + '\t')
#        f.write(str(e_out_array[i]) + '\t' + str(a_out_array_AU[i]) + '\t')
#        f.write(str(r1_array[i]) + '\t' + str(r2_array[i]) + '\t' + str(r3_array[i]) + '\t')
#        f.write('\n')
#    f.close()


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

    plt.savefig('plots/orbit/TPS_inner_orbit'+generic_name+'.pdf')
    plt.show()




    plt.plot(times_array_Myr,e_in_array)
    plt.plot(times_array_Myr,e_in_array, '.')
    plt.xlim(0,t_max_Myr)
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$e_\mathrm{in}$')
    plt.savefig('plots/orbit/e_in_time_'+generic_name+'.pdf')
    plt.show()

#    a_in_final_theory =  a_in_array_AU[0] * (m1_array[0] * m2_array[0] / m1_array / m2_array)**2 #stable mt
    a_in_final_theory =  a_in_array_AU[0] * (m1_array[0] + m2_array[0]) / (m1_array + m2_array)#wind 
    plt.plot(times_array_Myr,a_in_array_AU)
    plt.plot(times_array_Myr,a_in_array_AU, '.')
    plt.plot(times_array_Myr,a_in_final_theory)
    plt.plot(times_array_Myr,a_in_final_theory, '.')
    
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$a_\mathrm{in}$')
    plt.savefig('plots/orbit/semi_in_time_'+generic_name+'.pdf')
    plt.show()
#    exit(0)


    constants = 6283728.92847 # constants.G #in au and solar mass per Myr
    corot_spin_inner = 1./np.sqrt(a_in_array_AU**3/(m1_array+m2_array))*constants
    corot_spin_outer = 1./np.sqrt(a_out_array_AU**3/(m1_array+m2_array+m3_array))*constants

    plt.plot(times_array_Myr,spin1_array, 'b-')
    plt.plot(times_array_Myr,spin1_array, 'b.')
    plt.plot(times_array_Myr,spin2_array, 'g-')
    plt.plot(times_array_Myr,spin2_array, 'g.')
    plt.plot(times_array_Myr,spin3_array, 'r-')
    plt.plot(times_array_Myr,spin3_array, 'r.')

    plt.plot(times_array_Myr,corot_spin_inner, 'c-')
    plt.plot(times_array_Myr,corot_spin_inner, 'c,')
    plt.plot(times_array_Myr,corot_spin_outer, 'm-')
    plt.plot(times_array_Myr,corot_spin_outer, 'm,')
    
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$spin$')
    plt.savefig('plots/orbit/spin_time_'+generic_name+'.pdf')
    plt.show()


    J_orb2 = m1_array**2 * m2_array**2 / (m1_array+m2_array) * a_in_array_AU * ( 1-e_in_array**2) #*G
    J_orb = np.sqrt(J_orb2)
    J_spin1 =  spin1_array * m1_array* r1_array**2 #gyration radius  
    J_spin2 =  spin2_array * m2_array* r2_array**2 #gyration radius  
    J_spin3 =  spin3_array * m3_array* r3_array**2 #  gyration radius

    plt.plot(times_array_Myr, J_orb)
    plt.plot(times_array_Myr,J_orb, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$J orbit$')
    plt.savefig('plots/orbit/Jorbit_time_'+generic_name+'.pdf')
    plt.show()

    plt.plot(times_array_Myr,J_spin1)
    plt.plot(times_array_Myr,J_spin1, '.')
    plt.plot(times_array_Myr,J_spin2)
    plt.plot(times_array_Myr,J_spin2, '.')
    plt.plot(times_array_Myr,J_spin3)
    plt.plot(times_array_Myr,J_spin3, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$J spin$')
    plt.savefig('plots/orbit/Jspin_time_'+generic_name+'.pdf')
    plt.show()
    
    plt.semilogy(times_array_Myr,gy1_array)
    plt.semilogy(times_array_Myr,gy1_array, '.')
    plt.semilogy(times_array_Myr,gy2_array)
    plt.semilogy(times_array_Myr,gy2_array, '.')
    plt.semilogy(times_array_Myr,gy3_array)
    plt.semilogy(times_array_Myr,gy3_array, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$gyration radius$')
    plt.savefig('plots/orbit/gyration_radius_time_'+generic_name+'.pdf')
    plt.show()
      
    

    plt.semilogy(times_array_Myr,r1_array)
    plt.semilogy(times_array_Myr,r1_array, '.')
    plt.semilogy(times_array_Myr,r2_array)
    plt.semilogy(times_array_Myr,r2_array, '.')
    plt.semilogy(times_array_Myr,r3_array)
    plt.semilogy(times_array_Myr,r3_array, '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$radius$')
    plt.savefig('plots/orbit/radius_time_'+generic_name+'.pdf')
    plt.show()
    
#    dr1_array =r1_array[1:]-r1_array[:-1]
#    dr2_array =r2_array[1:]-r2_array[:-1]
#    dr3_array =r3_array[1:]-r3_array[:-1]
#    dt_array = times_array_Myr[1:] - times_array_Myr[:-1]
#    plt.semilogy(times_array_Myr[1:], dr1_array/dt_array)
#    plt.semilogy(times_array_Myr[1:], dr1_array/dt_array, '.')
#    plt.semilogy(times_array_Myr[1:], dr2_array/dt_array)
#    plt.semilogy(times_array_Myr[1:], dr2_array/dt_array, '.')
#    plt.semilogy(times_array_Myr[1:], dr3_array/dt_array)
#    plt.semilogy(times_array_Myr[1:], dr3_array/dt_array, '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$dr/dt$')
#    plt.savefig('plots/orbit/drdt_time_'+generic_name+'.pdf')
#    plt.show()
#
#
#
#
#    plt.semilogy(times_array_Myr[1:], dr1_array/r1_array[1:])
#    plt.semilogy(times_array_Myr[1:], dr1_array/r1_array[1:], '.')
#    plt.semilogy(times_array_Myr[1:], dr2_array/r2_array[1:])
#    plt.semilogy(times_array_Myr[1:], dr2_array/r2_array[1:], '.')
#    plt.semilogy(times_array_Myr[1:], dr3_array/r3_array[1:])
#    plt.semilogy(times_array_Myr[1:], dr3_array/r3_array[1:], '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$dr/r$')
#    plt.savefig('plots/orbit/dr_time_'+generic_name+'.pdf')
#    plt.show()



    dm1_array =m1_array[:-1]-m1_array[1:]
    dm2_array =m2_array[:-1]-m2_array[1:]
    dm3_array =m3_array[:-1]-m3_array[1:]
    
    plt.plot(times_array_Myr[1:], m1_array[1:])
    plt.plot(times_array_Myr[1:], m1_array[1:], '.')
    plt.plot(times_array_Myr[1:], m2_array[1:])
    plt.plot(times_array_Myr[1:], m2_array[1:], '.')
    plt.plot(times_array_Myr[1:], m3_array[1:])
    plt.plot(times_array_Myr[1:], m3_array[1:], '.')
    plt.xlabel('$t/\mathrm{Myr}$')
    plt.ylabel('$mass$')
    plt.savefig('plots/orbit/mass_time_'+generic_name+'.pdf')
    plt.show()

#
#    plt.semilogy(times_array_Myr[1:], dm1_array/dt_array)
#    plt.semilogy(times_array_Myr[1:], dm1_array/dt_array, '.')
#    plt.semilogy(times_array_Myr[1:], dm2_array/dt_array)
#    plt.semilogy(times_array_Myr[1:], dm2_array/dt_array, '.')
#    plt.semilogy(times_array_Myr[1:], dm3_array/dt_array)
#    plt.semilogy(times_array_Myr[1:], dm3_array/dt_array, '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$dm/dt$')
#    plt.savefig('plots/orbit/dmdt_time_'+generic_name+'.pdf')
#    plt.show()
#
#    plt.semilogy(times_array_Myr[1:], dm1_array/m1_array[1:])
#    plt.semilogy(times_array_Myr[1:], dm1_array/m1_array[1:], '.')
#    plt.semilogy(times_array_Myr[1:], dm2_array/m2_array[1:])
#    plt.semilogy(times_array_Myr[1:], dm2_array/m2_array[1:], '.')
#    plt.semilogy(times_array_Myr[1:], dm3_array/m3_array[1:])
#    plt.semilogy(times_array_Myr[1:], dm3_array/m3_array[1:], '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$dm/m$')
#    plt.savefig('plots/orbit/dm_time_'+generic_name+'.pdf')
#    plt.show()
#
# 
#    plt.semilogy(times_array_Myr[1:], dt_array)
#    plt.semilogy(times_array_Myr[1:], dt_array, '.')
#    plt.semilogy(times_array_Myr[1:], dt_array)
#    plt.semilogy(times_array_Myr[1:], dt_array, '.')
#    plt.semilogy(times_array_Myr[1:], dt_array)
#    plt.semilogy(times_array_Myr[1:], dt_array, '.')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$dt$')
#    plt.savefig('plots/orbit/dt_time_'+generic_name+'.pdf')
#    plt.show()
#




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
#    plt.savefig('plots/orbit/semi_inner_wind'+generic_name+'.pdf')
#    plt.show()
#    
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot/a_in_array_AU)
#    plt.plot(times_array_Myr, a_in_array_AU[0]*Mtot[0]/Mtot/a_in_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{in}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig('plots/orbit/semi_inner_rel_wind'+generic_name+'.pdf')
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
#    plt.savefig('plots/orbit/semi_inner_cons_mt'+generic_name+'.pdf')
#    plt.show()
#    
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2/a_in_array_AU)
#    plt.plot(times_array_Myr, a_in_array_AU[0]*(m1_array[0]*m2_array[0]/m1_array/m2_array)**2/a_in_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{in}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig('plots/orbit/semi_inner_rel_cons_mt'+generic_name+'.pdf')
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
    plt.savefig('plots/orbit/TPS_outer_orbit'+generic_name+'.pdf')
    plt.show()



#    plt.plot(times_array_Myr,e_out_array)
#    plt.xlim(0,t_max_Myr)
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.ylabel('$e_\mathrm{out}$')
#    plt.savefig('plots/orbit/e_out_time_'+generic_name+'.pdf')
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
#    plt.savefig('plots/orbit/semi_outer_wind'+generic_name+'.pdf')
#    plt.show()
#
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott/a_out_array_AU)
#    plt.plot(times_array_Myr, a_out_array_AU[0]*Mtott[0]/Mtott/a_out_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{out}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig('plots/orbit/semi_outer_rel_wind'+generic_name+'.pdf')
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
#    plt.savefig('plots/orbit/semi_outer_cons_mt'+generic_name+'.pdf')
#    plt.show()
#
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2/a_out_array_AU)
#    plt.plot(times_array_Myr, a_out_array_AU[0]*(m_in_array[0]*m3_array[0]/m_in_array/m3_array)**2/a_out_array_AU, '.')
#    plt.ylabel('$relative error a_\mathrm{out}$')
#    plt.xlabel('$t/\mathrm{Myr}$')
#    plt.savefig('plots/orbit/semi_outer_rel_cons_mt'+generic_name+'.pdf')
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
            tend= 5.0 |units.Myr, number = 0, maximum_radius_change_factor = 0.005,
            tidal_terms = True):


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
            metallicity, tend, number, maximum_radius_change_factor, tidal_terms)


    triple_class_object.evolve_model()
#    plot_function(triple_class_object)
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
    parser.add_option("-r", dest="maximum_radius_change_factor", type="float", default = 0.005,
                      help="maximum_radius_change_factor [%default] %unit")
    parser.add_option("--tidal", dest="tidal_terms", action="store_false",
                      help="tidal terms included [%default] %unit")
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
    plot_function(triple_class_object)

    if REPORT_TRIPLE_EVOLUTION:
        print 'succesfully finished'


