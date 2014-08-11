## TPS:         Triple population synthesis.
##              computes the evolution of a population of given triples   
##              given any initial conditions (M, q, m, A, a, E, e, i, G, g, O, o).
## Output:      in the form of the following files:
##              not yet implemented
## Options:   --M_max    upper limit for the inner primary mass [100 Msun]
##            --M_min    lower limit for the inner primary mass [0.1 Msun]
##            --M_distr  mass function option: 0) Kroupa [default]
##                                             1) Scalo
##                                             2) Miller & Scalo
##                                             3) Salpeter
##                                             4) Logarithmically flat 
##                                             5) Eggleton
##            --Q_max    upper limit for the inner mass ratio [1.]
##            --Q_min    lower limit for the inner mass ratio [0.]
##            --Q_distr  inner mass ratio option: 0) Flat (uniform) distribution [default]
##                                                1) Kroupa IMF
##            --q_max    upper limit for the outer mass ratio [1.]
##            --q_min    lower limit for the mass of the outer star [0.]
##            --q_distr  outer mass ratio option: 0) Flat (uniform) distribution [default]
##                                                1) Kroupa IMF
##            --A_max    upper limit for the inner semi-major axis [5e6 RSun]
##            --A_min    lower limit for the inner semi-major axis [5]
##            --A_distr  inner semi-major axis option: 0) logFlat distribution [default]
##                                                     1) Constant semi-major axis
##            --a_max    upper limit for the outer semi-major axis [5e6 RSun]
##            --a_min    lower limit for the outer semi-major axis [5 RSun]
##            --a_distr  outer semi-major axis option: 0) logFlat distribution [default]
##                                                     1) Constant semi-major axis
##            --E_max    upper limit for the inner eccentricity [1.]
##            --E_min    lower limit for the inner eccentricity [0.]
##            --E_distr  inner eccentricity option: 0) Thermal [default]
##                                                  1) Constant eccentricity
##            --e_max    upper limit for the outer eccentricity [1.]
##            --e_min    lower limit for the outer eccentricity [0.]
##            --e_distr  outer eccentricity option: 0) Thermal [default]
##                                                  1) Constant eccentricity
##            --i_max    upper limit for the mutual inclination [pi]
##            --i_min    lower limit for the mutual inclination [0]
##            --i_distr  mutual inclination option: 0) Circular uniform distribution [default]
##                                                  1) Constant inclination
##            --G_max    upper limit for the inner argument of pericenter []
##            --G_min    lower limit for the inner argument of pericenter []
##            --G_distr  inner argument of pericenter option: 0) Uniform [default]
##                                                            1) Constant inclination
##            --g_max    upper limit for the outer argument of pericenter []
##            --g_min    lower limit for the outer argument of pericenter []
##            --g_distr  outer argument of pericenter option: 0) Uniform [default]
##                                                            1) Constant inclination
##            --O_max    upper limit for the inner longitude of ascending node []
##            --O_min    lower limit for the inner longitude of ascending node []
##            --O_distr  inner longitude of ascending node option: 0) ? [default]
##            --o_max    upper limit for the outer longitude of ascending node []
##            --o_min    lower limit for the outer longitude of ascending node []
##            --o_distr  outer longitude of ascending node option: 0) ? [default]
##            -T or -t   binary end time. [13500 Myr]
##            -z         metallicity of stars  [0.02 Solar] 
##            -n         number of binaries to be simulated.  [1]

#not implemented yet
##            -D        stopping condition at merger or disruption [False]
##            -s         random seed


import triple

from amuse.units.optparse import OptionParser
from amuse.units import units, constants
from amuse.support.console import set_printing_strategy
import numpy as np

from amuse.ic.kroupa import new_kroupa_mass_distribution
from amuse.ic.scalo import new_scalo_mass_distribution
from amuse.ic.millerscalo import new_miller_scalo_mass_distribution
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.flatimf import new_flat_mass_distribution


min_mass = 0.1 |units.MSun # for stars
max_mass = 100 |units.MSun



def flat_distr(lower, upper):
    if lower.unit != upper.unit:
        print 'flat_distr: different units'
        exit(1)    
    return np.random.uniform(lower.number, upper.number)|lower.unit

def log_flat_distr(lower, upper):
    if lower.unit != upper.unit:
        print 'log_flat_distr: different units'
        exit(1)    
    return 10**np.random.uniform(np.log10(lower.number), np.log10(upper.number))|lower.unit

def thermal_distr(lower, upper): #unit is missing for eccentricity
#    if lower.unit != upper.unit:
#        print 'thermal_distr: different units'
#        exit(1)    
#    return np.sqrt(np.random.uniform(lower.number, upper.number))|lower.unit
    return np.sqrt(np.random.uniform(lower, upper))

def circular_uniform_distr(lower, upper): #unit is missing for inclination
#    print lower, upper
#    if lower.unit != upper.unit:
#        print 'circular_uniform_distr: different units'
#        exit(1)    
#    return np.arccos(np.random.uniform(np.cos(lower.number), np.cos(upper.number)))|lower.unit

    return np.arccos(np.random.uniform(np.cos(lower), np.cos(upper)))

def eggleton_mass_distr(lower_mass, upper_mass):
    turnover_mass = 0.3|units.MSun
    power = 0.85
    
    y_max = (upper_mass/turnover_mass) **(1/power)
    upper = y_max / (1+y_max)
    y_min = (lower_mass/turnover_mass) **(1/power)
    lower = y_min / (1+y_min)

    x = np.random.uniform(lower, upper)
    y=turnover_mass * (x/(1-x))**power
    return turnover_mass * (x/(1-x))**power


class Generate_initial_triple:
    #-------
    #setup stellar system
    def __init__(self, in_primary_mass_max, in_primary_mass_min, 
                        in_mass_ratio_max, in_mass_ratio_min, 
                        out_mass_ratio_max, out_mass_ratio_min, 
                        in_semi_max, in_semi_min, out_semi_max, out_semi_min, 
                        in_ecc_max, in_ecc_min, out_ecc_max, out_ecc_min, 
                        incl_max, incl_min,
                        in_aop_max, in_aop_min, out_aop_max, out_aop_min,
                        in_loan_max, in_loan_min, out_loan_max, out_loan_min,
                        in_primary_mass_distr, in_mass_ratio_distr, out_mass_ratio_distr,
                        in_semi_distr,  out_semi_distr, in_ecc_distr, out_ecc_distr, incl_distr,
                        in_aop_distr, out_aop_distr, in_loan_distr, out_loan_distr):
                            
                        self.generate_mass(in_primary_mass_max, in_primary_mass_min, 
                            in_mass_ratio_max, in_mass_ratio_min,
                            out_mass_ratio_max, out_mass_ratio_min,
                            in_primary_mass_distr, in_mass_ratio_distr, out_mass_ratio_distr)
                            
                        self.generate_semi(in_semi_max, in_semi_min, 
                            out_semi_max, out_semi_min,
                            in_semi_distr,  out_semi_distr)

                        self.generate_ecc(in_ecc_max, in_ecc_min, 
                            out_ecc_max, out_ecc_min,
                            in_ecc_distr, out_ecc_distr)

                        self.generate_incl(incl_max, incl_min, incl_distr)

                        self.generate_aop(in_aop_max, in_aop_min, 
                            out_aop_max, out_aop_min,
                            in_aop_distr, out_aop_distr)

                        self.generate_loan(in_loan_max, in_loan_min,                   
                            out_loan_max, out_loan_min,
                            in_loan_distr, out_loan_distr)

    #-------                        
    def generate_mass(self, in_primary_mass_max, in_primary_mass_min, 
                        in_mass_ratio_max, in_mass_ratio_min,
                        out_mass_ratio_max, out_mass_ratio_min,
                        in_primary_mass_distr, in_mass_ratio_distr, out_mass_ratio_distr):

        if in_primary_mass_max == in_primary_mass_min:
            self.in_primary_mass = in_primary_mass_min
        else:
            if in_primary_mass_distr == 1: #Scalo 1986
                self.in_primary_mass= new_scalo_mass_distribution(1, in_primary_mass_max)[0]
                while self.in_primary_mass < in_primary_mass_min:
                        self.in_primary_mass = new_scalo_mass_distribution(1, in_primary_mass_max)[0]
            elif in_primary_mass_distr == 2:#Miller & Scale 1979
                self.in_primary_mass= new_miller_scalo_mass_distribution(1, in_primary_mass_max)[0]
                while self.in_primary_mass < in_primary_mass_min:
                        self.in_primary_mass = new_miller_scalo_mass_distribution(1, in_primary_mass_max)[0]
            elif in_primary_mass_distr == 3: #Salpeter with slope 2.35
                self.in_primary_mass= new_salpeter_mass_distribution(1, in_primary_mass_min, in_primary_mass_max)[0]
            elif in_primary_mass_distr == 4: # Flat in log space
                self.in_primary_mass= new_flat_mass_distribution(1, in_primary_mass_min, in_primary_mass_max)[0]
            elif in_primary_mass_distr == 5: # Eggleton 2009, 399, 1471, Salpeter-like with turnover at low masses
                self.in_primary_mass= eggleton_mass_distr(in_primary_mass_min, in_primary_mass_max)
            else: #Kroupa 2001
                self.in_primary_mass = new_kroupa_mass_distribution(1, in_primary_mass_max)[0]
                while self.in_primary_mass < in_primary_mass_min:
                        self.in_primary_mass = new_kroupa_mass_distribution(1, in_primary_mass_max)[0]


        if in_mass_ratio_max == in_mass_ratio_min:
            in_mass_ratio = in_mass_ratio_min
            self.in_secondary_mass = in_mass_ratio * self.in_primary_mass
        else: 
            if in_mass_ratio_distr == 1:# Kroupa 2001 
                self.in_secondary_mass = new_kroupa_mass_distribution(1, in_primary_mass_max)[0]
                while self.in_secondary_mass < in_primary_mass_min:
                        self.in_secondary_mass = new_kroupa_mass_distribution(1, in_primary_mass_max)[0]
            else: # flat distribution
               in_mass_ratio = flat_distr(in_mass_ratio_min, in_mass_ratio_max)
               self.in_secondary_mass = in_mass_ratio * self.in_primary_mass        


        if out_mass_ratio_max == out_mass_ratio_min:
            out_mass_ratio = out_mass_ratio_min
            self.out_mass = out_mass_ratio * (self.in_primary_mass + self.in_secondary_mass)
        else: 
            if out_mass_ratio_distr == 1:# Kroupa 2001 
                self.out_mass = new_kroupa_mass_distribution(1, in_primary_mass_max)[0]
                while self.out_mass < in_primary_mass_min:
                        self.out_mass = new_kroupa_mass_distribution(1, in_primary_mass_max)[0]
            else: # flat distribution
               out_mass_ratio = flat_distr(out_mass_ratio_min, out_mass_ratio_max)
               self.out_mass = out_mass_ratio * (self.in_primary_mass + self.in_secondary_mass)
               print out_mass_ratio, self.out_mass


    def generate_semi(self, 
                        in_semi_max, in_semi_min, 
                        out_semi_max, out_semi_min,
                        in_semi_distr,  out_semi_distr):
                        
        if in_semi_max == in_semi_min:
            self.in_semi = in_semi_min
        else:
            if in_semi_distr == 1: #Constant 
                print 'TPS::generate_semi: unambiguous choise of constant semi-major axis'
                print '--A_min option to set the value of the semi-major axis in the inner binary'                
                self.in_ecc = in_semi_min
            else: # log flat distribution
                 maximal_semi = min(in_semi_max, out_semi_max)
                 self.in_semi = log_flat_distr(in_semi_min, maximal_semi)
                        
                        
        if out_semi_max == out_semi_min:
            self.out_semi = out_semi_min
        else:
            if out_semi_distr == 1: #Constant 
                print 'TPS::generate_semi: unambiguous choise of constant semi-major axis'
                print '--a_min option to set the value of the semi-major axis in the outer binary'                
                self.out_ecc = out_semi_min
            else: # log flat distribution
                 minimal_semi = max(out_semi_min, self.in_semi) # outer orbit is always larger then inner orbit
                 self.out_semi = log_flat_distr(minimal_semi, out_semi_max)
                         
 
        if self.out_semi < self.in_semi:
            print self.in_semi, self.out_semi
            print in_semi_min, in_semi_max
            print out_semi_min, out_semi_max
            print 'error generate_semi: a_out< a_in'
            exit(1)

    def generate_ecc(self,
                        in_ecc_max, in_ecc_min, 
                        out_ecc_max, out_ecc_min,
                        in_ecc_distr, out_ecc_distr):
                        
        if in_ecc_max == in_ecc_min:
            self.in_ecc = in_ecc_min
        else:
            if in_ecc_distr == 1: #Constant 
                print 'TPS::generate_ecc: unambiguous choise of constant eccentricity'
                print '--E_min option to set the value of the eccentricity in the inner binary'                
                self.in_ecc = in_ecc_min
            else: #Thermal distribution
                 self.in_ecc = thermal_distr(in_ecc_min, in_ecc_max)
                 
                 
        if out_ecc_max == out_ecc_min:
            self.out_ecc = out_ecc_min
        else:
            if out_ecc_distr == 1: #Constant 
                print 'TPS::generate_ecc: unambiguous choise of constant eccentricity'
                print '--e_min option to set the value of eccentricity in the outer binary'                
                self.out_ecc = out_ecc_min
            else: #Thermal distribution
                 self.out_ecc = thermal_distr(out_ecc_min, out_ecc_max)


    def generate_incl(self, incl_max, incl_min, incl_distr):
        if incl_max == incl_min:
            self.incl = incl_min
        else:
            if incl_distr == 1: #Constant 
                print 'TPS::generate_incl: unambiguous choise of constant mutual inclination'
                print '--i_min option to set the value of the mutual inclination in the inner triple'                
                self.incl = incl_min
            else: #Circular uniform distribution
                 self.incl = circular_uniform_distr(incl_min, incl_max)

    def generate_aop(self,
                        in_aop_max, in_aop_min, 
                        out_aop_max, out_aop_min,
                        in_aop_distr, out_aop_distr):

        if in_aop_max == in_aop_min:
            self.in_aop = in_aop_min
        else:
            if in_aop_distr == 1: #Constant 
                print 'TPS::generate_aop: unambiguous choise of constant argument of pericenter'
                print '--G_min option to set the value of the argument of pericenter of the inner binary'                
                self.in_aop = in_aop_min
            else: #Uniform distribution
#                 self.in_aop = flat_distr(in_aop_min, in_aop_max)
                 self.in_aop = np.random.uniform(in_aop_min, in_aop_max) # unit is missing
                 
        if out_aop_max == out_aop_min:
            self.out_aop = out_aop_min
        else:
            if out_aop_distr == 1: #Constant 
                print 'TPS::generate_aop: unambiguous choise of constant argument of pericenter'
                print '--g_min option to set the value of the argument of pericenter of the outer binary'                
                self.out_aop = out_aop_min
            else: #Uniform distribution
#                 self.out_aop = flat_distr(out_aop_min, out_aop_max)
                 self.out_aop = np.random.uniform(out_aop_min, out_aop_max) # unit is missing
                 

    def generate_loan(self,
                        in_loan_max, in_loan_min,                   
                        out_loan_max, out_loan_min,
                        in_loan_distr, out_loan_distr):                                
        self.in_loan = 0.
        self.out_loan = 0.
#-------
        
#-------
    def print_triple(self):
        print '\nTriple - ' 
        print 'm =', self.in_primary_mass, self.in_secondary_mass, self.out_mass
        print 'a =', self.in_semi, self.out_semi
        print 'e =', self.in_ecc, self.out_ecc
        print 'i =', self.incl
        print 'g =', self.in_aop, self.out_aop
        print 'o =', self.in_loan, self.out_loan 
#-------

#-------
def evolve_triples(in_primary_mass_max, in_primary_mass_min, 
                        in_mass_ratio_max, in_mass_ratio_min, 
                        out_mass_ratio_max, out_mass_ratio_min, 
                        in_semi_max, in_semi_min, out_semi_max, out_semi_min, 
                        in_ecc_max, in_ecc_min, out_ecc_max, out_ecc_min, 
                        incl_max, incl_min,
                        in_aop_max, in_aop_min, out_aop_max, out_aop_min,
                        in_loan_max, in_loan_min, out_loan_max, out_loan_min, 
                        in_primary_mass_distr, in_mass_ratio_distr, out_mass_ratio_distr,
                        in_semi_distr,  out_semi_distr, in_ecc_distr, out_ecc_distr, incl_distr,
                        in_aop_distr, out_aop_distr, in_loan_distr, out_loan_distr,                                                                     
                        metallicity, tend, number, stop_at_merger_or_disruption, seed):


    for i_n in range(number):
        triple_system = Generate_initial_triple(in_primary_mass_max, in_primary_mass_min, 
                    in_mass_ratio_max, in_mass_ratio_min, 
                    out_mass_ratio_max, out_mass_ratio_min, 
                    in_semi_max, in_semi_min, out_semi_max, out_semi_min, 
                    in_ecc_max, in_ecc_min, out_ecc_max, out_ecc_min, 
                    incl_max, incl_min,
                    in_aop_max, in_aop_min, out_aop_max, out_aop_min,
                    in_loan_max, in_loan_min, out_loan_max, out_loan_min,
                    in_primary_mass_distr, in_mass_ratio_distr, out_mass_ratio_distr,
                    in_semi_distr,  out_semi_distr, in_ecc_distr, out_ecc_distr, incl_distr,
                    in_aop_distr, out_aop_distr, in_loan_distr, out_loan_distr)
        
        triple_system.print_triple()
        
#        triple.main(inner_primary_mass = triple_system.in_primary_mass, 
#                    inner_secondary_mass = triple_system.in_secondary_mass, 
#                    outer_mass = triple_system.out_mass, 
#                    inner_semimajor_axis = triple_system.in_semi, 
#                    outer_semimajor_axis = triple_system.out_semi, 
#                    inner_eccentricity = triple_system.in_ecc, 
#                    outer_eccentricity = triple_system.out_ecc, 
#                    mutual_inclination = triple_system.incl, 
#                    inner_argument_of_pericenter = triple_system.in_aop, 
#                    outer_argument_of_pericenter = triple_system.out_aop, 
#                    inner_longitude_of_ascending_node = triple_system.in_loan, 
#                    outer_longitude_of_ascending_node = triple_system.out_loan, 
#                    metallicity = metallicity, tend = tend)                        
                                

def test_initial_parameters(in_primary_mass_max, in_primary_mass_min, 
                        in_mass_ratio_max, in_mass_ratio_min, 
                        out_mass_ratio_max, out_mass_ratio_min, 
                        in_semi_max, in_semi_min, out_semi_max, out_semi_min, 
                        in_ecc_max, in_ecc_min,  out_ecc_max, out_ecc_min, 
                        incl_max, incl_min,
                        in_aop_max, in_aop_min, out_aop_max, out_aop_min,
                        in_loan_max, in_loan_min, out_loan_max, out_loan_min,   
                        in_primary_mass_distr, in_mass_ratio_distr, out_mass_ratio_distr,
                        in_semi_distr,  out_semi_distr, in_ecc_distr, out_ecc_distr, incl_distr,
                        in_aop_distr, out_aop_distr, in_loan_distr, out_loan_distr,                   
                        metallicity, tend, number, stop_at_merger_or_disruption, seed):

    if (in_primary_mass_min < min_mass) or (in_primary_mass_max > max_mass):
        print 'error: inner primary mass not in allowed range [', min_mass, ',', max_mass, ']'
        exit(1)
        
    if (in_primary_mass_max < in_primary_mass_min):
        print 'error: maximum inner primary mass smaller than minimum in primary mass'
        exit(1)

        
    if (in_mass_ratio_min < 0.) or (in_mass_ratio_max > 1.):
        print 'error: inner mass ratio not in allowed range'
        exit(1)
        
    if (in_mass_ratio_max < in_mass_ratio_min):
        print 'error: maximum inner mass ratio smaller than minimum mass ratio'
        exit(1)


    if (out_mass_ratio_min < 0.) or (out_mass_ratio_max > 1.):
        print 'error: outer mass ratio not in allowed range'
        exit(1)
        
    if (out_mass_ratio_max < out_mass_ratio_min):
        print 'error: maximum outer mass ratio smaller than minimum mass ratio'
        exit(1)
        

    if (in_semi_min < 5.|units.RSun) or (in_semi_max > 5.e6|units.RSun):
        print 'error: inner separation not in allowed range [5,5e6]RSun'
        exit(1)

    if (out_semi_min < 5.|units.RSun) or (in_semi_max > 5.e6|units.RSun):
        print 'error: outer separation not in allowed range [5,5e6]RSun'
        exit(1)
        
    if (in_semi_max < in_semi_min):
        print 'error: maximum inner separation smaller than minimum in separation'
        exit(1)

    if (out_semi_max < out_semi_min):
        print 'error: maximum outer separation smaller than minimum outer separation'
        exit(1)
        
    if (in_semi_min > out_semi_max):
        print 'error: maximum outer separation smaller than minimum inner separation - no overlap for inner and outer orbit'
        exit(1)
        
        


    if (in_ecc_min < 0.) or (in_ecc_max > 1.):
        print 'error: inner eccentricity not in allowed range [0,1]'
        exit(1)

    if (out_ecc_min < 0.) or (out_ecc_max > 1.):
        print 'error: outer eccentricity not in allowed range [0,1]'
        exit(1)

    if (in_ecc_max < in_ecc_min):
        print 'error: maximum inner eccentricity smaller than minimum ecc'
        exit(1)

    if (out_ecc_max < out_ecc_min):
        print 'error: maximum outer eccentricity smaller than minimum ecc'
        exit(1)


    if (incl_min < 0.) or (incl_max > np.pi):
        print 'error: mutual inclination not in allowed range [0, pi]'
        exit(1)

    if (incl_max < incl_min):
        print 'error: maximum mutual inclination smaller than minimum mutual inclination'
        exit(1)



    if (in_aop_min < 0.) or (in_aop_max > 2*np.pi):
        print 'error: inner argument of pericenter not in allowed range [0,2*pi]'
        exit(1)

    if (out_aop_min < 0.) or (out_aop_max > 2*np.pi):
        print 'error: outer argument of pericenter not in allowed range [0,2*pi]'
        exit(1)

    if (in_aop_max < in_aop_min):
        print 'error: maximum inner argument of pericenter smaller than minimum argument of pericenter'
        exit(1)

    if (out_aop_max < out_aop_min):
        print 'error: maximum outer argument of pericenter smaller than minimum argument of pericenter'
        exit(1)


    if (in_loan_min < 0.) or (in_loan_max > 2*np.pi):
        print 'error: inner longitude of ascending node not in allowed range [0,2*pi]'
        exit(1)

    if (out_loan_min < 0.) or (out_loan_max > 2*np.pi):
        print 'error: outer longitude of ascending node not in allowed range [0,2*pi]'
        exit(1)

    if (in_loan_max < in_loan_min):
        print 'error: maximum inner longitude of ascending node smaller than minimum argument of pericenter'
        exit(1)

    if (out_loan_max < out_loan_min):
        print 'error: maximum outer longitude of ascending node smaller than minimum argument of pericenter'
        exit(1)



def parse_arguments():
    parser = OptionParser()
    parser.add_option("--M_min", unit=units.MSun, 
                      dest="in_primary_mass_min", type="float", default = 1.|units.MSun,
                      help="minimum of inner primary mass [%default]")
    parser.add_option("--M_max", unit=units.MSun, 
                      dest="in_primary_mass_max", type="float", default = 100|units.MSun,
                      help="maximum of inner primary mass [%default]")
    parser.add_option("--M_distr", dest="in_primary_mass_distr", type="int", default = 0,
                      help="inner primary mass distribution [Kroupa]")
                      
    parser.add_option("--Q_max", unit=units.none, 
                      dest="in_mass_ratio_max", type="float", default = 1.0 |units.none,
                      help="maximum of inner mass ratio [%default]")
    parser.add_option("--Q_min", unit=units.none, 
                      dest="in_mass_ratio_min", type="float", default = 0. |units.none,
                      help="minimum of inner mass ratio [%default]")
    parser.add_option("--Q_distr", dest="in_mass_ratio_distr", type="int", default = 0,
                      help="inner mass ratio distribution [Flat]")

    parser.add_option("--q_max", unit=units.none, 
                      dest="out_mass_ratio_max", type="float", default = 1.0 |units.none,
                      help="maximum of outer mass ratio [%default]")
    parser.add_option("--q_min", unit=units.none, 
                      dest="out_mass_ratio_min", type="float", default = 0. |units.none,
                      help="minimum of outer mass ratio [%default]")
    parser.add_option("--q_distr", dest="out_mass_ratio_distr", type="int", default = 0,
                      help="outer mass ratio distribution [Flat]")
                      

    parser.add_option("--A_min", unit=units.RSun,
                      dest="in_semi_min", type="float", 
                      default = 5|units.RSun,
                      help="minimum of inner semi major axis [%default]")
    parser.add_option("--A_max", unit=units.RSun,
                      dest="in_semi_max", type="float", 
                      default = 5e6|units.RSun,
                      help="maximum of inner semi major axis [%default]")
    parser.add_option("--A_distr", dest="in_semi_distr", type="int", default = 0,
                      help="inner semimajor axis distribution [?]")


    parser.add_option("--a_min", unit=units.RSun,
                      dest="out_semi_min", type="float", 
                      default = 5|units.RSun,
                      help="minimum of outer semi major axis [%default]")
    parser.add_option("--a_max", unit=units.RSun,
                      dest="out_semi_max", type="float", 
                      default = 5e6|units.RSun,
                      help="maximum of outer semi major axis [%default]")
    parser.add_option("--a_distr", dest="out_semi_distr", type="int", default = 0,
                      help="outer semimajor axis distribution [?]")

    parser.add_option("--E_min",
                      dest="in_ecc_min", type="float", default = 0.,
                      help="minimum of inner eccentricity [%default]")
    parser.add_option("--E_max",
                      dest="in_ecc_max", type="float", default = 1.0,
                      help="maximum of inner eccentricity [%default]")
    parser.add_option("--E_distr", dest="in_ecc_distr", type="int", default = 0,
                      help="inner eccentricity distribution [Thermal]")

    parser.add_option("--e_min",
                      dest="out_ecc_min", type="float", default = 0.,
                      help="minimum of outer eccentricity [%default]")
    parser.add_option("--e_max",
                      dest="out_ecc_max", type="float", default = 1.,
                      help="maximum of outer eccentricity [%default]")
    parser.add_option("--e_distr", dest="out_ecc_distr", type="int", default = 0,
                      help="outer eccentricity distribution [Thermal]")
                      
                      
    parser.add_option("--i_min, --I_min",
                      dest="incl_min", type="float", default = 0.0,
                      help="minimum of mutual inclination [rad] [%default]")
    parser.add_option("--i_max, --I_max",
                      dest="incl_max", type="float", default = np.pi,
                      help="maximum of mutual inclination [rad] [%default]")
    parser.add_option("--i_distr, --I_distr", dest="incl_distr", type="int", default = 0,
                      help="mutual inclination distribution [Circular uniform]")

                      
    parser.add_option("--G_min",
                      dest="in_aop_min", type="float", default = 0.,
                      help="minimum of inner argument of pericenter [rad] [%default]")
    parser.add_option("--G_max",
                      dest="in_aop_max", type="float", default = 2*np.pi,
                      help="maximum of inner argument of pericenter [rad] [%default]")
    parser.add_option("--G_distr", dest="in_aop_distr", type="int", default = 0,
                      help="inner argument of pericenter distribution [Uniform]")

    parser.add_option("--g_min",
                      dest="out_aop_min", type="float", default = 0.5,
                      help="minimum of outer argument of pericenter [rad] [%default]")
    parser.add_option("--g_max",
                      dest="out_aop_max", type="float", default = 0.5,
                      help="maximum of outer argument of pericenter [rad] [%default]")
    parser.add_option("--g_distr", dest="out_aop_distr", type="int", default = 0,
                      help="outer argument of pericenter distribution [Uniform]")
                      
    parser.add_option("--O_min",
                      dest="in_loan_min", type="float", default = 0,
                      help="minimum of inner longitude of ascending node [rad] [%default]")
    parser.add_option("--O_max",
                      dest="in_loan_max", type="float", default = 0,
                      help="maximum of inner longitude of ascending node [rad] [%default]")
    parser.add_option("--O_distr", dest="in_loan_distr", type="int", default = 0,
                      help="inner longitude of ascending node distribution [?]")

    parser.add_option("--o_min",
                      dest="out_loan_min", type="float", default = 0,
                      help="minimum of outer longitude of ascending node [rad] [%default]")
    parser.add_option("--o_max",
                      dest="out_loan_max", type="float", default = 0,
                      help="maximum of outer longitude of ascending node [rad] [%default]")
    parser.add_option("--o_distr", dest="out_loan_distr", type="int", default = 0,
                      help="outer longitude of ascending node distribution [?]")
                      
    parser.add_option("-t", "-T", unit=units.Myr, 
                      dest="tend", type="float", default = 200|units.Myr,
                      help="end time [%default] %unit")
    parser.add_option("-z", unit=units.none, 
                      dest="metallicity", type="float", default = 0.02|units.none,
                      help="metallicity [%default] %unit")
    parser.add_option("-n", dest="number", type="int", default = 10,
                      help="number of systems [%default]")
    parser.add_option("-D", unit=units.none, 
                      dest="stop_at_merger_or_disruption", type="int", default = 0|units.none,
                      help="stop at merger or disruption [%default] %unit") #should be bool
    parser.add_option("-s", unit=units.none, 
                      dest="seed", type="float", default = 0.|units.none,
                      help="seed [%default] %unit")
#    int actual_seed = srandinter(input_seed);



    options, args = parser.parse_args()
    return options.__dict__



if __name__ == '__main__':
    options = parse_arguments()
    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    test_initial_parameters(**options)
    evolve_triples(**options)



