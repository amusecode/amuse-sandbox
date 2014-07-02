from amuse.lab import *
from binary import *
from math import sqrt

### added by Adrian ###
from amuse.community.seculartriple_TPS.interface import SecularTriple
import numpy as np
from amuse.units import quantities
from matplotlib import pyplot
#######################

#constants
timestep_factor = 0.01
#error_dm
#error_dr

REPORT_TRIPLE_EVOLUTION = False

class Triple:
    #-------
    def __init__(self, inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            mutual_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node,
            metallicity,
            tend):

        stars = self.make_stars(inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis)
        bins = self.make_bins(stars, inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node)

        triples = Particles(1)
        triples[0].inner_binary = bins[0]
        triples[0].outer_binary = bins[1]
        triples[0].mutual_inclination = mutual_inclination

        self.is_star = False
        self.is_binary = False
        self.is_triple = True
        self.first_contact = True
        self.tend = tend #...
        self.time = 0.0|units.yr
        self.previous_time = 0.0|units.yr

        self.particles = triples
        self.setup_se_code(metallicity, stars)
        self.setup_secular_code(triples)

        self.update_previous_se_parameters()
        self.update_se_parameters() 

    def make_stars(self, inner_primary_mass, inner_secondary_mass, outer_mass, inner_semimajor_axis, outer_semimajor_axis):
        stars = Particles(3)
        stars.is_star = True
        stars.is_binary = False 
        stars.is_donor = False

        stars[0].mass = inner_primary_mass
        stars[1].mass = inner_secondary_mass
        stars[2].mass = outer_mass

        # default now corotating
        corotating_angular_frequency_inner = 1./np.sqrt(inner_semimajor_axis**3/constants.G / (stars[0].mass + stars[1].mass))
        corotating_angular_frequency_outer = 1./np.sqrt(outer_semimajor_axis**3/constants.G / (stars[0].mass + stars[1].mass + stars[2].mass))
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
        
        bins[0].child1 = stars[0]
        bins[0].child2 = stars[1]
        bins[0].semimajor_axis = inner_semimajor_axis
        bins[0].eccentricity = inner_eccentricity
        bins[0].argument_of_pericenter = inner_argument_of_pericenter
        bins[0].longitude_of_ascending_node = inner_longitude_of_ascending_node
        
        bins[0].mass_transfer_rate = 0.0 | units.MSun/units.yr
        bins[0].accretion_efficiency_mass_transfer = 1.0

        bins[1].semimajor_axis = outer_semimajor_axis
        bins[1].eccentricity = outer_eccentricity
        bins[1].child1 = stars[2]
        bins[1].child2 = bins[0]
        bins[1].argument_of_pericenter = outer_argument_of_pericenter                
        bins[1].longitude_of_ascending_node = outer_longitude_of_ascending_node
        
        bins[1].mass_transfer_rate = 0.0 | units.MSun/units.yr        
        bins[1].accretion_efficiency_mass_transfer = 1.0

        # binary evolutionary settings
        bins[0].specific_AM_loss_mass_transfer = 2.5 
        bins[1].specific_AM_loss_mass_transfer = 2.5

        
        for x in bins:
            x.mass = x.child1.mass + x.child2.mass

        return bins
    #-------
            
    #-------
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
        self.secular_code.parameters.include_inner_tidal_terms = False
        self.secular_code.parameters.include_outer_tidal_terms = False
        self.secular_code.parameters.include_inner_wind_terms = False
        self.secular_code.parameters.include_outer_wind_terms = False
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



    def determine_timestep(self):         
        if REPORT_TRIPLE_EVOLUTION:
            print "Dt=", self.se_code.particles.time_step, self.tend/100.0
    
        ### for testing/plotting purposes only ###
        timestep = self.tend/100.0 
    
        timestep = min(timestep, 
            self.particles[0].outer_binary.child1.time_step,
            self.particles[0].inner_binary.child1.time_step,
            self.particles[0].inner_binary.child2.time_step)
            
        #during stable mass transfer     
        if (self.particles[0].inner_binary.child1.is_donor):
            timestep = min(timestep, abs(timestep_factor*self.particles[0].inner_binary.child1.mass/self.particles[0].inner_binary.child1.mass_transfer_rate))
        if (self.particles[0].inner_binary.child2.is_donor):
            timestep = min(timestep, abs(timestep_factor*self.particles[0].inner_binary.child2.mass/self.particles[0].inner_binary.child2.mass_transfer_rate))
        if (self.particles[0].outer_binary.child1.is_donor):
            timestep = min(timestep, abs(timestep_factor*self.particles[0].outer_binary.child1.mass/self.particles[0].outer_binary.child1.mass_transfer_rate))
            
        self.time += timestep
    
    #-------
    def update_previous_se_parameters(self):
        self.previous_time = self.time

        self.particles[0].inner_binary.child1.previous_mass = self.particles[0].inner_binary.child1.mass 
        self.particles[0].inner_binary.child2.previous_mass = self.particles[0].inner_binary.child2.mass 
        self.particles[0].outer_binary.child1.previous_mass = self.particles[0].outer_binary.child1.mass 
        self.particles[0].outer_binary.child2.previous_mass = self.particles[0].outer_binary.child2.mass 

        self.particles[0].inner_binary.child1.previous_radius = self.particles[0].inner_binary.child1.radius 
        self.particles[0].inner_binary.child2.previous_radius = self.particles[0].inner_binary.child2.radius 
        self.particles[0].outer_binary.child1.previous_radius = self.particles[0].outer_binary.child1.radius 
    #-------

    #-------
    def update_se_wind_parameters(self):
        self.update_wind_mass_loss_rate()
        self.update_time_derivative_of_radius()
        self.update_envelope_mass()
                    
    def update_wind_mass_loss_rate(self):
        #update wind mass loss rate
        #note wind mass loss rate < 0
        timestep = self.time - self.previous_time
        if timestep > 0|units.yr: 
            self.particles[0].inner_binary.child1.wind_mass_loss_rate = (self.particles[0].inner_binary.child1.mass - self.particles[0].inner_binary.child1.previous_mass)/timestep
            self.particles[0].inner_binary.child2.wind_mass_loss_rate = (self.particles[0].inner_binary.child2.mass - self.particles[0].inner_binary.child2.previous_mass)/timestep
            self.particles[0].outer_binary.child1.wind_mass_loss_rate = (self.particles[0].outer_binary.child1.mass - self.particles[0].outer_binary.child1.previous_mass)/timestep
        else:
            #initialization
            self.particles[0].inner_binary.child1.wind_mass_loss_rate = 0.0|units.MSun/units.yr
            self.particles[0].inner_binary.child2.wind_mass_loss_rate = 0.0|units.MSun/units.yr
            self.particles[0].outer_binary.child1.wind_mass_loss_rate = 0.0|units.MSun/units.yr
            
    def update_time_derivative_of_radius(self):
        #update time_derivative_of_radius for effect of wind on spin
        #radius change due to stellar evolution, not mass transfer
        timestep = self.time - self.previous_time
        if timestep > 0|units.yr:
            self.particles[0].inner_binary.child1.time_derivative_of_radius = (self.particles[0].inner_binary.child1.radius - self.particles[0].inner_binary.child1.previous_radius)/timestep
            self.particles[0].inner_binary.child2.time_derivative_of_radius = (self.particles[0].inner_binary.child2.radius - self.particles[0].inner_binary.child2.previous_radius)/timestep
            self.particles[0].outer_binary.child1.time_derivative_of_radius = (self.particles[0].outer_binary.child1.radius - self.particles[0].outer_binary.child1.previous_radius)/timestep
        else:
            #initialization
            self.particles[0].inner_binary.child1.time_derivative_of_radius = 0.0 | units.RSun/units.yr
            self.particles[0].inner_binary.child2.time_derivative_of_radius = 0.0 | units.RSun/units.yr
            self.particles[0].outer_binary.child1.time_derivative_of_radius = 0.0 | units.RSun/units.yr

    #-------

    #-------
    def update_se_parameters(self):
        self.update_convective_envelope_mass()
        self.update_convective_envelope_radius()
        self.update_gyration_radius()
        

    def update_convective_envelope_mass(self):
        #update convective envelope radius
        #the prescription of Hurley, Pols & Tout 2000 is implemented in SeBa, however note that the prescription in BSE is different
        self.particles[0].inner_binary.child1.convective_envelope_mass = self.particles[0].inner_binary.child1.get_convective_envelope_mass()
        self.particles[0].inner_binary.child2.convective_envelope_mass = self.particles[0].inner_binary.child2.get_convective_envelope_mass()
        self.particles[0].outer_binary.child1.convective_envelope_mass = self.particles[0].outer_binary.child1.get_convective_envelope_mass()

        
    def update_convective_envelope_radius(self):
        #update convective envelope radius
        #the prescription of Hurley, Tout & Pols 2002 is implemented in SeBa, , however note that the prescription in BSE is different 
        self.particles[0].inner_binary.child1.convective_envelope_radius = self.particles[0].inner_binary.child1.get_convective_envelope_radius() 
        self.particles[0].inner_binary.child2.convective_envelope_radius = self.particles[0].inner_binary.child2.get_convective_envelope_radius() 
        self.particles[0].outer_binary.child1.convective_envelope_radius = self.particles[0].outer_binary.child1.get_convective_envelope_radius()

    def update_gyration_radius(self):
        #update gyration radius
        self.particles[0].inner_binary.child1.gyration_radius = self.se_code.particles[0].get_gyration_radius_sq()**0.5
        self.particles[0].inner_binary.child2.gyration_radius = self.se_code.particles[1].get_gyration_radius_sq()**0.5
        self.particles[0].outer_binary.child1.gyration_radius = self.se_code.particles[2].get_gyration_radius_sq()**0.5
    #-------
        
        

        

    
def evolve_center_of_mass(binary):
    print "evolve center of mass"
    exit(-1)
    return


def resolve_triple_interaction(triple):

    if REPORT_TRIPLE_EVOLUTION:
        print '\ninner binary'
    if triple.particles[0].inner_binary.is_binary:
        resolve_binary_interaction(triple.particles[0].inner_binary, triple)
    elif triple.particles[0].inner_binary.is_star:
#        'e.g. if system merged'
        evolve_center_of_mass(triple.particles[0].inner_binary, triple)
    else:
        print 'resolve triple interaction: type of inner system unknown'
        exit(-1)                    

    if REPORT_TRIPLE_EVOLUTION:
        print '\nouter binary'
    if triple.particles[0].outer_binary.is_binary:
        resolve_binary_interaction(triple.particles[0].outer_binary, triple)
    elif triple.particles[0].outer_binary.is_star:
#        'e.g. if there is no outer star?'
        evolve_center_of_mass(triple.particles[0].outer_binary, triple)
    else:
        print 'resolve triple interaction: type of outer system unknown'
        exit(-1)                    



def safety_check_timestep(triple):
        dm_1 = (triple.particles[0].inner_binary.child1.mass - triple.particles[0].inner_binary.child1.previous_mass)/triple.particles[0].inner_binary.child1.mass
        dm_2 = (triple.particles[0].inner_binary.child2.mass - triple.particles[0].inner_binary.child2.previous_mass)/triple.particles[0].inner_binary.child2.mass
        dm_3 = (triple.particles[0].outer_binary.child1.mass - triple.particles[0].outer_binary.child1.previous_mass)/triple.particles[0].outer_binary.child1.mass
                
        if REPORT_TRIPLE_EVOLUTION:
            print 'wind mass loss rates:', 
            print triple.particles[0].inner_binary.child1.wind_mass_loss_rate,
            print triple.particles[0].inner_binary.child2.wind_mass_loss_rate,
            print triple.particles[0].outer_binary.child1.wind_mass_loss_rate
            print 'relative wind mass losses:',
            print dm_1, dm_2, dm_3
        error_dm = 0.05
        if (dm_1 > error_dm) or (dm_2 > error_dm) or (dm_3 > error_dm):
            print 'Change in mass in a single timestep larger then', error_dm
            print dm_1, dm_2, dm_3
            print triple.particles[0].inner_binary.child1.stellar_type
            print triple.particles[0].inner_binary.child2.stellar_type
            print triple.particles[0].outer_binary.child1.stellar_type
            exit(-1)


        dr_1 = (triple.particles[0].inner_binary.child1.radius - triple.particles[0].inner_binary.child1.previous_radius)/triple.particles[0].inner_binary.child1.radius
        dr_2 = (triple.particles[0].inner_binary.child2.radius - triple.particles[0].inner_binary.child2.previous_radius)/triple.particles[0].inner_binary.child2.radius
        dr_3 = (triple.particles[0].outer_binary.child1.radius - triple.particles[0].outer_binary.child1.previous_radius)/triple.particles[0].outer_binary.child1.radius
        if REPORT_TRIPLE_EVOLUTION:    
            print 'change in radius over time:', 
            print triple.particles[0].inner_binary.child1.time_derivative_of_radius,
            print triple.particles[0].inner_binary.child2.time_derivative_of_radius,
            print triple.particles[0].outer_binary.child1.time_derivative_of_radius
            print 'relative change in radius:',
            print dr_1, dr_2, dr_3
        error_dr = 0.05
        if (dr_1 > error_dr) or (dr_2 > error_dr) or (dr_3 > error_dr):
            print 'Change in radius in a single timestep larger then', error_dr
            print dr_1, dr_2, dr_3
            print triple.particles[0].inner_binary.child1.stellar_type
            print triple.particles[0].inner_binary.child2.stellar_type
            print triple.particles[0].outer_binary.child1.stellar_type
            exit(-1)
        

    
def evolve_triple(triple):
    
    # for plotting data
    times_array = quantities.AdaptingVectorQuantity() 
    a_in_array = quantities.AdaptingVectorQuantity()
    e_in_array = []
    i_mutual_array = []
    g_in_array = []


    while triple.time<triple.tend:
        triple.update_previous_se_parameters()
        triple.determine_timestep()        

        # update stellar parameters in the secular code
        triple.channel_to_secular.copy()   
        
        # do secular evolution
        if triple.is_triple == True:
            triple.secular_code.evolve_model(triple.time)
        else:
            print 'Secular code disabled'
            exit(-1)

        if triple.secular_code.triples[0].dynamical_instability == True:
            print "Dynamical instability at time/Myr = ",time.value_in(units.Myr)
            exit(0)
        if triple.secular_code.triples[0].inner_collision == True:
            print "Inner collision at time/Myr = ",time.value_in(units.Myr)
            exit(0)
        if triple.secular_code.triples[0].outer_collision == True:
            print "Outer collision at time/Myr = ",time.value_in(units.Myr)
            exit(0)
        # this updates orbital elements within triple.particles[0] 
        triple.channel_from_secular.copy() 


        # when the secular code finds that mass transfer starts, evolve the stars only until that time
        if ((triple.particles[0].inner_binary.child1.is_donor or
            triple.particles[0].inner_binary.child2.is_donor or
            triple.particles[0].outer_binary.child1.is_donor) and
            triple.first_contact):
                if REPORT_TRIPLE_EVOLUTION:
                    print 'Times:', triple.previous_time, triple.time, triple.secular_code.model_time
                triple.time = triple.secular_code.model_time 
                
                # mass should not be transferred just yet -> check if this works ok
                # do not overwrite this parameter in the secular code    
                triple.particles[0].inner_binary.child1.is_donor = False
                triple.particles[0].inner_binary.child2.is_donor = False
                triple.particles[0].outer_binary.child1.is_donor = False


        #do stellar evolution 
        triple.se_code.evolve_model(triple.time)
        triple.channel_from_se.copy()
        triple.update_se_wind_parameters()
        
        safety_check_timestep(triple)
        
        ##what should be the order, first mass transfer or first stellar evolution?
        resolve_triple_interaction(triple)        
#        Rl1, Rl2, Rl3 = triple.secular_code.give_roche_radii(triple.particles[0])
        triple.update_se_parameters()
        
        #should also do safety check timestep here
        
        
        # for plotting data
        times_array.append(triple.time)
        e_in_array.append(triple.particles[0].inner_binary.eccentricity)
        a_in_array.append(triple.particles[0].inner_binary.semimajor_axis)
        g_in_array.append(triple.particles[0].inner_binary.argument_of_pericenter)    
        i_mutual_array.append(triple.particles[0].mutual_inclination)        
        
    # for plotting data
    e_in_array = np.array(e_in_array)
    g_in_array = np.array(g_in_array)
    i_mutual_array = np.array(i_mutual_array)
    triple.plot_data = plot_data_container()
    triple.plot_data.times_array = times_array
    triple.plot_data.a_in_array = a_in_array
    triple.plot_data.e_in_array = e_in_array
    triple.plot_data.i_mutual_array = i_mutual_array
    triple.plot_data.g_in_array = g_in_array
    #########################################
    
class plot_data_container():
    def __init__(self):
        return

def plot_function(triple):
    ### plots to test secular code ###
    ##################################
    
    figure = pyplot.figure(figsize=(10,13))
    N_subplots = 4

    plot_e_in = figure.add_subplot(N_subplots,1,1)
    plot_i_mutual = figure.add_subplot(N_subplots,1,2)
    plot_e_in_g_in = figure.add_subplot(N_subplots,1,3)
    plot_a_in = figure.add_subplot(N_subplots,1,4)

    times_array_Myr = triple.plot_data.times_array.value_in(units.Myr)
    a_in_array_AU = triple.plot_data.a_in_array.value_in(units.AU)
    g_in_array = triple.plot_data.g_in_array
    e_in_array = triple.plot_data.e_in_array
    i_mutual_array = triple.plot_data.i_mutual_array
    plot_e_in.plot(times_array_Myr,e_in_array)
    plot_i_mutual.plot(times_array_Myr,i_mutual_array*180.0/np.pi)
    plot_e_in_g_in.plot(np.cos(g_in_array),e_in_array)
    plot_a_in.plot(times_array_Myr,a_in_array_AU)

    t_max_Myr = max(times_array_Myr)
    plot_e_in.set_xlim(0,t_max_Myr)
    plot_i_mutual.set_xlim(0,t_max_Myr)
    plot_i_mutual.set_ylim(0.9*min(i_mutual_array*180.0/np.pi),1.1*max(i_mutual_array*180.0/np.pi))

    plot_e_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_mutual.set_xlabel('$t/\mathrm{Myr}$')
    plot_e_in_g_in.set_xlabel('$\cos(g_\mathrm{in})$')

    plot_e_in.set_ylabel('$e_\mathrm{in}$')
    plot_i_mutual.set_ylabel('$i_\mathrm{mutual} ({}^\circ)$')
    plot_e_in_g_in.set_ylabel('$e_\mathrm{in}$')
    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)

    pyplot.show()




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
                      dest="inner_longitude_of_ascending_node", type="float", default = 0,
                      help="inner longitude of ascending node [rad] [%default]")
    parser.add_option("-O",
                      dest="outer_longitude_of_ascending_node", type="float", default = 0,
                      help="outer longitude of ascending node [rad] [%default]")
    parser.add_option("-M", unit=units.MSun, 
                      dest="inner_primary_mass", type="float", default = 1.3|units.MSun,
                      help="inner primary mass [%default]")
    parser.add_option("-m",  unit=units.MSun, 
                      dest="inner_secondary_mass", type="float", default = 0.5|units.MSun,
                      help="inner secondary mass [%default]")
    parser.add_option("-n",  unit=units.MSun, 
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
    evolve_triple(triple)

    plot_function(triple)
