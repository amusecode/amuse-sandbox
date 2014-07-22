import triple

from amuse.units.optparse import OptionParser
from amuse.units import units, constants
from amuse.support.console import set_printing_strategy
import numpy as np



def initialize_triples(inner_primary_mass_max, inner_primary_mass_min, 
                        inner_secondary_mass_max, inner_secondary_mass_min, 
                        outer_mass_max, outer_mass_min,
                        inner_semimajor_axis_max, inner_semimajor_axis_min, 
                        outer_semimajor_axis_max, outer_semimajor_axis_min, 
                        inner_eccentricity_max, inner_eccentricity_min, 
                        outer_eccentricity_max, outer_eccentricity_min, 
                        mutual_inclination_max, mutual_inclination_min,
                        inner_argument_of_pericenter_max, inner_argument_of_pericenter_min, 
                        outer_argument_of_pericenter_max, outer_argument_of_pericenter_min,
                        inner_longitude_of_ascending_node_max, inner_longitude_of_ascending_node_min,         
                        outer_longitude_of_ascending_node_max, outer_longitude_of_ascending_node_min,                        
                        metallicity, tend, 
                        number, stop_at_merger_or_disruption):
    
    print 'no populations yet'
    print number
    print 'grid based or monte carlo based' 
    
#    tr = triple.main(inner_primary_mass = inner_primary_mass, inner_secondary_mass = inner_secondary_mass, outer_mass = outer_mass, inner_semimajor_axis = inner_semimajor_axis_max, outer_semimajor_axis = outer_semimajor_axis_max, inner_eccentricity = inner_eccentricity_min, outer_eccentricity = outer_eccentricity_min, mutual_inclination = mutual_inclination_min, inner_argument_of_pericenter = inner_argument_of_pericenter_min, outer_argument_of_pericenter = outer_argument_of_pericenter_min, inner_longitude_of_ascending_node = inner_longitude_of_ascending_node_min, outer_longitude_of_ascending_node = outer_longitude_of_ascending_node_min, tend = tend)

    print inner_semimajor_axis_min


def parse_arguments():
    parser = OptionParser()
    parser.add_option("--a_min", unit=units.AU,
                      dest="inner_semimajor_axis_min", type="float", 
                      default = 100|units.AU,
                      help="minimum of inner semi major axis [%default]")
    parser.add_option("--a_max", unit=units.AU,
                      dest="inner_semimajor_axis_max", type="float", 
                      default = 100|units.AU,
                      help="maximum of inner semi major axis [%default]")

    parser.add_option("--A_min", unit=units.AU,
                      dest="outer_semimajor_axis_min", type="float", 
                      default = 10000|units.AU,
                      help="minimum of outer semi major axis [%default]")
    parser.add_option("--A_max", unit=units.AU,
                      dest="outer_semimajor_axis_max", type="float", 
                      default = 10000|units.AU,
                      help="maximum of outer semi major axis [%default]")

    parser.add_option("--e_min",
                      dest="inner_eccentricity_min", type="float", default = 0,
                      help="minimum of inner eccentricity [%default]")
    parser.add_option("--e_max",
                      dest="inner_eccentricity_max", type="float", default = 0,
                      help="maximum of inner eccentricity [%default]")

    parser.add_option("--E_min",
                      dest="outer_eccentricity_min", type="float", default = 0,
                      help="minimum of outer eccentricity [%default]")
    parser.add_option("--E_max",
                      dest="outer_eccentricity_max", type="float", default = 0,
                      help="maximum of outer eccentricity [%default]")
                      
                      
    parser.add_option("--M_min", unit=units.MSun, 
                      dest="inner_primary_mass_min", type="float", default = 2|units.MSun,
                      help="minimum of inner primary mass [%default]")
    parser.add_option("--M_max", unit=units.MSun, 
                      dest="inner_primary_mass_max", type="float", default = 2|units.MSun,
                      help="maximum of inner primary mass [%default]")
    parser.add_option("--m_min",  unit=units.MSun, 
                      dest="inner_secondary_mass_min", type="float", default = 1|units.MSun,
                      help="minimum of inner secondary mass [%default]")
    parser.add_option("--m_max",  unit=units.MSun, 
                      dest="inner_secondary_mass_max", type="float", default = 1|units.MSun,
                      help="maximum of inner secondary mass [%default]")
    parser.add_option("--l_min",  unit=units.MSun, 
                      dest="outer_mass_min", type="float", default = 0.5|units.MSun,
                      help="minimum of outer mass [%default]")
    parser.add_option("--l_max",  unit=units.MSun, 
                      dest="outer_mass_max", type="float", default = 0.5|units.MSun,
                      help="maximum of outer mass [%default]")
                      
    parser.add_option("--i_min",
                      dest="mutual_inclination_min", type="float", default = 80.0*np.pi/180.0,
                      help="minimum of mutual inclination [rad] [%default]")
    parser.add_option("--i_max",
                      dest="mutual_inclination_max", type="float", default = 80.0*np.pi/180.0,
                      help="maximum of mutual inclination [rad] [%default]")
                      
    parser.add_option("--g_min",
                      dest="inner_argument_of_pericenter_min", type="float", default = 0.1,
                      help="minimum of inner argument of pericenter [rad] [%default]")
    parser.add_option("--g_max",
                      dest="inner_argument_of_pericenter_max", type="float", default = 0.1,
                      help="maximum of inner argument of pericenter [rad] [%default]")

    parser.add_option("--G_min",
                      dest="outer_argument_of_pericenter_min", type="float", default = 0.5,
                      help="minimum of outer argument of pericenter [rad] [%default]")
    parser.add_option("--G_max",
                      dest="outer_argument_of_pericenter_max", type="float", default = 0.5,
                      help="maximum of outer argument of pericenter [rad] [%default]")
                      
    parser.add_option("--o_min",
                      dest="inner_longitude_of_ascending_node_min", type="float", default = 0,
                      help="minimum of inner longitude of ascending node [rad] [%default]")
    parser.add_option("--o_max",
                      dest="inner_longitude_of_ascending_node_max", type="float", default = 0,
                      help="maximum of inner longitude of ascending node [rad] [%default]")

    parser.add_option("--O_min",
                      dest="outer_longitude_of_ascending_node_min", type="float", default = 0,
                      help="minimum of outer longitude of ascending node [rad] [%default]")
    parser.add_option("--O_max",
                      dest="outer_longitude_of_ascending_node_max", type="float", default = 0,
                      help="maximum of outer longitude of ascending node [rad] [%default]")
                      
    parser.add_option("-t", "-T", unit=units.Myr, 
                      dest="tend", type="float", default = 200|units.Myr,
                      help="end time [%default] %unit")
    parser.add_option("-z", unit=units.none, 
                      dest="metallicity", type="float", default = 0.02|units.none,
                      help="metallicity [%default] %unit")

    parser.add_option("-n",  unit=units.none, 
                      dest="number", type="int", default = 10|units.none,
                      help="number of systems [%default]")
    parser.add_option("-D", unit=units.none, 
                      dest="stop_at_merger_or_disruption", type="int", default = 0|units.none,
                      help="stop at merger or disruption [%default] %unit") #should be bool



    options, args = parser.parse_args()
    return options.__dict__


if __name__ == '__main__':
    options = parse_arguments()
    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    initialize_triples(**options)



