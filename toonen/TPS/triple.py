from amuse.lab import *
from binary import *

### added by Adrian ###
from amuse.community.seculartriple_TPS.interface import SecularTriple
import numpy
from amuse.units import quantities
from matplotlib import pyplot
#######################

REPORT_TRIPLE_EVOLUTION = False

class Triple:
    def __init__(self, inner_primary_mass, inner_secondary_mass, outer_mass,
            inner_semimajor_axis, outer_semimajor_axis,
            inner_eccentricity, outer_eccentricity,
            mutual_inclination,
            inner_argument_of_pericenter, outer_argument_of_pericenter,
            inner_longitude_of_ascending_node, outer_longitude_of_ascending_node,
            metallicity,
            tend):

        stars = self.make_stars(inner_primary_mass, inner_secondary_mass, outer_mass)
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
        self.tend = tend #...
        self.time = 0.0|units.yr
        self.previous_time = 0.0|units.yr

        self.particles = triples
        self.setup_se_code(metallicity, stars)
        self.setup_secular_code(triples)

#        self.update_previous()
#        self.update() 

    def make_stars(self, inner_primary_mass, inner_secondary_mass, outer_mass):
        stars = Particles(3)
        stars.is_star = True
        stars.is_binary = False 
        stars.is_donor = False

        stars[0].mass = inner_primary_mass
        stars[1].mass = inner_secondary_mass
        stars[2].mass = outer_mass

        stars[0].wind_mass_loss_rate = -1.0e-6 | units.MSun/units.yr
        stars[1].wind_mass_loss_rate = 0.0 | units.MSun/units.yr
        stars[2].wind_mass_loss_rate = 0.0 | units.MSun/units.yr

        stars[0].spin_angular_frequency = 0.0 | 1.0/units.yr
        stars[1].spin_angular_frequency = 0.0 | 1.0/units.yr
        stars[2].spin_angular_frequency = 0.0 | 1.0/units.yr

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
        bins[0].accretion_efficiency_wind_child1_to_child2 = 0.0
        bins[0].accretion_efficiency_wind_child2_to_child1 = 0.0
        bins[0].accretion_efficiency_mass_transfer = 1.0
        bins[0].specific_AM_loss_mass_transfer = 2.5

        bins[1].semimajor_axis = outer_semimajor_axis
        bins[1].eccentricity = outer_eccentricity
        bins[1].child1 = stars[2]
        bins[1].child2 = bins[0]
        bins[1].argument_of_pericenter = outer_argument_of_pericenter                
        bins[1].longitude_of_ascending_node = outer_longitude_of_ascending_node
        
        bins[1].mass_transfer_rate = 0.0 | units.MSun/units.yr        
        bins[1].accretion_efficiency_wind_child1_to_child2 = 0.0
        bins[1].accretion_efficiency_wind_child2_to_child1 = 0.0
        bins[1].accretion_efficiency_mass_transfer = 1.0
        bins[1].specific_AM_loss_mass_transfer = 2.5
        
        for x in bins:
            x.mass = x.child1.mass + x.child2.mass

        return bins
        
    def setup_se_code(self, metallicity, stars):
        self.se_code = SeBa()
        self.se_code.parameters.metallicity = metallicity
        self.se_code.particles.add_particles(stars)
        self.channel_from_se = self.se_code.particles.new_channel_to(stars)
        self.channel_to_se = stars.new_channel_to(self.se_code.particles)
        self.channel_from_se.copy_attributes(["mass", "radius", "stellar_type", "core_mass"])
    
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

    def update_previous(self):
        self.previous_time = self.time

        self.particles[0].inner_binary.child1.previous_mass = self.particles[0].inner_binary.child1.mass 
        self.particles[0].inner_binary.child2.previous_mass = self.particles[0].inner_binary.child2.mass 
        self.particles[0].outer_binary.child1.previous_mass = self.particles[0].outer_binary.child1.mass 
        self.particles[0].outer_binary.child2.previous_mass = self.particles[0].outer_binary.child2.mass 


    def update(self):
        #update envelope mass
        self.particles[0].inner_binary.child1.envelope_mass = self.particles[0].inner_binary.child1.mass - self.particles[0].inner_binary.child1.core_mass
        self.particles[0].inner_binary.child2.envelope_mass = self.particles[0].inner_binary.child2.mass - self.particles[0].inner_binary.child2.core_mass
        self.particles[0].outer_binary.child1.envelope_mass = self.particles[0].outer_binary.child1.mass - self.particles[0].outer_binary.child1.core_mass
        
        #update wind mass loss rate
        timestep = self.previous_time - self.time
        if timestep > 0|units.yr: #maybe better to get the rate directly out of seba
            self.particles[0].inner_binary.child1.wind_mass_loss_rate = (self.particles[0].inner_binary.child1.previous_mass - self.particles[0].inner_binary.child1.mass)/timestep
            self.particles[0].inner_binary.child2.wind_mass_loss_rate = (self.particles[0].inner_binary.child2.mass - self.particles[0].inner_binary.child2.previous_mass)/timestep
            self.particles[0].outer_binary.child1.wind_mass_loss_rate = (self.particles[0].outer_binary.child1.mass - self.particles[0].outer_binary.child1.previous_mass)/timestep
        else:
            #initialization
            self.particles[0].inner_binary.child1.wind_mass_loss_rate = 0.0|units.MSun/units.yr
            self.particles[0].inner_binary.child2.wind_mass_loss_rate = 0.0|units.MSun/units.yr
            self.particles[0].outer_binary.child1.wind_mass_loss_rate = 0.0|units.MSun/units.yr

        #update gyration radius
        self.particles[0].inner_binary.child1.gyration_radius = self.se_code.particles[0].get_gyration_radius_sq()**0.5
        self.particles[0].inner_binary.child2.gyration_radius = self.se_code.particles[1].get_gyration_radius_sq()**0.5
        self.particles[0].outer_binary.child1.gyration_radius = self.se_code.particles[2].get_gyration_radius_sq()**0.5

    
def evolve_center_of_mass(binary):
    print "evolve center of mass"
    exit(-1)
    return


def resolve_triple_interaction(triple):

    if REPORT_TRIPLE_EVOLUTION:
        print '\ninner binary'
    if triple.particles[0].inner_binary.is_binary:
        resolve_binary_interaction(triple.particles[0].inner_binary)
    elif triple.particles[0].inner_binary.is_star:
#        'e.g. if system merged'
        evolve_center_of_mass(triple.particles[0].inner_binary)
    else:
        print 'resolve triple interaction: type of inner system unknown'
        exit(-1)                    

    if REPORT_TRIPLE_EVOLUTION:
        print '\nouter binary'
    if triple.particles[0].outer_binary.is_binary:
        resolve_binary_interaction(triple.particles[0].outer_binary)
    elif triple.particles[0].outer_binary.is_star:
#        'e.g. if there is no outer star?'
        evolve_center_of_mass(triple.particles[0].outer_binary)
    else:
        print 'resolve triple interaction: type of outer system unknown'
        exit(-1)                    

def evolve_triple(triple):
    
    ### for testing/plotting purposes only ###
    triple.timestep = triple.tend/100.0 
    ##########################################
    
    ### temporary; only for plotting data ###
    times_array = quantities.AdaptingVectorQuantity() 
    a_in_array = quantities.AdaptingVectorQuantity()
    e_in_array = []
    i_mutual_array = []
    g_in_array = []
    #########################################

    while triple.time<triple.tend:
        triple.time += triple.timestep
        
        if REPORT_TRIPLE_EVOLUTION:
            print "Dt=", triple.se_code.particles.age, triple.time, triple.se_code.particles.time_step

        triple.update_previous()
        #eventually adjustable timestep, minimum of stellar and secular evolution        
        triple.se_code.evolve_model(triple.time)
        triple.channel_from_se.copy()
        triple.update()


        ### do secular triple evolution ###
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
        triple.channel_from_secular.copy() ### this updates orbital elements within triple.particles[0] ###
        ###################################

        resolve_triple_interaction(triple)        
#        triple.channel_to_se.copy()#masses
        
        ### temporary; only for plotting data ###
        times_array.append(triple.time)
        e_in_array.append(triple.particles[0].inner_binary.eccentricity)
        a_in_array.append(triple.particles[0].inner_binary.semimajor_axis)
        g_in_array.append(triple.particles[0].inner_binary.argument_of_pericenter)    
        i_mutual_array.append(triple.particles[0].mutual_inclination)        
        #########################################
        
    ### temporary; only for plotting data ###
    e_in_array = numpy.array(e_in_array)
    g_in_array = numpy.array(g_in_array)
    i_mutual_array = numpy.array(i_mutual_array)
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
    plot_i_mutual.plot(times_array_Myr,i_mutual_array*180.0/numpy.pi)
    plot_e_in_g_in.plot(numpy.cos(g_in_array),e_in_array)
    plot_a_in.plot(times_array_Myr,a_in_array_AU)

    t_max_Myr = max(times_array_Myr)
    plot_e_in.set_xlim(0,t_max_Myr)
    plot_i_mutual.set_xlim(0,t_max_Myr)
    plot_i_mutual.set_ylim(0.9*min(i_mutual_array*180.0/numpy.pi),1.1*max(i_mutual_array*180.0/numpy.pi))

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
                      dest="mutual_inclination", type="float", default = 80.0*numpy.pi/180.0,
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
