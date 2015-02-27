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
##                                                     2) Tokovinin
##            --a_max    upper limit for the outer semi-major axis [5e6 RSun]
##            --a_min    lower limit for the outer semi-major axis [5 RSun]
##            --a_distr  outer semi-major axis option: 0) logFlat distribution [default]
##                                                     1) Constant semi-major axis
##                                                     2) Tokovinin
##            --E_max    upper limit for the inner eccentricity [1.]
##            --E_min    lower limit for the inner eccentricity [0.]
##            --E_distr  inner eccentricity option: 0) Thermal [default]
##                                                  1) Constant eccentricity
##            --e_max    upper limit for the outer eccentricity [1.]
##            --e_min    lower limit for the outer eccentricity [0.]
##            --e_distr  outer eccentricity option: 0) Thermal [default]
##                                                  1) Constant eccentricity
##            --i_max    upper limit for the relative inclination [pi]
##            --i_min    lower limit for the relative inclination [0]
##            --i_distr  relative inclination option: 0) Circular uniform distribution [default]
##                                                  1) Constant inclination
##            --G_max    upper limit for the inner argument of pericenter [2pi]
##            --G_min    lower limit for the inner argument of pericenter [0]
##            --G_distr  inner argument of pericenter option: 0) Uniform distribution [default]
##                                                            1) Constant argument of pericenter
##            --g_max    upper limit for the outer argument of pericenter [2pi]
##            --g_min    lower limit for the outer argument of pericenter [0]
##            --g_distr  outer argument of pericenter option: 0) Uniform distribution [default]
##                                                            1) Constant argument of pericenter
##             outer longitude of ascending nodes = inner - pi               
##            --O_max    upper limit for the inner longitude of ascending node [pi]
##            --O_min    lower limit for the inner longitude of ascending node [0]
##            --O_distr  inner longitude of ascending node option: 0) Circular uniform distribution
##                                                            1) Constant longitude of ascending nodes [default]
##            -T or -t   binary end time. [13500 Myr]
##            -z         metallicity of stars  [0.02 Solar] 
##            -n         number of triples to be simulated.  [1]
##            -N         number of initial triple.  [0]
##            --stop_at_merger                  stopping condition at merger [True]
##            --stop_at_disintegrated           stopping condition at disintegration [True]
##            --stop_at_triple_mass_transfer    stopping condition at mass transfer in outer binary [True]
##            --stop_at_collision               stopping condition at collision [True]
##            --stop_at_dynamical_instability   stopping condition at dynamical instability [True]
##            --stop_at_mass_transfer           stopping condition at mass transfer [False]
         
#not implemented yet
##            -s         random seed


#import triple
import triple_nokozaidt as triple

from amuse.units.optparse import OptionParser
from amuse.units import units, constants
from amuse.support.console import set_printing_strategy
import numpy as np
from scipy.interpolate import interp1d

from amuse.ic.kroupa import new_kroupa_mass_distribution
from amuse.ic.scalo import new_scalo_mass_distribution
from amuse.ic.millerscalo import new_miller_scalo_mass_distribution
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.flatimf import new_flat_mass_distribution


min_mass = 0.08 |units.MSun # for stars
max_mass = 100 |units.MSun
REPORT = False 
REPORT_USER_WARNINGS = False 
stop_at_init_mass_transfer = True

def flat_distr(lower, upper):
    return np.random.uniform(lower, upper)

def log_flat_distr(lower, upper):
    lower_RSun = lower.value_in(units.RSun)
    upper_RSun = upper.value_in(units.RSun)
    x= np.random.uniform(np.log10(lower_RSun), np.log10(upper_RSun))
    return (10**x)|units.RSun 
    
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
                        in_loan_max, in_loan_min, 
                        in_primary_mass_distr, in_mass_ratio_distr, out_mass_ratio_distr,
                        in_semi_distr,  out_semi_distr, in_ecc_distr, out_ecc_distr, incl_distr,
                        in_aop_distr, out_aop_distr, in_loan_distr):
                        
                        if in_primary_mass_distr == 5:
                            convergence = False
                            while convergence == False:
                                convergence = self.generate_mass_and_semi_eggleton(in_primary_mass_max, in_primary_mass_min, in_semi_max, in_semi_min, 
                    out_semi_max, out_semi_min)                            
#                                   print convergence
                        else:    
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

                        self.generate_loan(in_loan_max, in_loan_min, in_loan_distr)

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
#                while self.in_secondary_mass < in_primary_mass_min:
#                        self.in_secondary_mass = new_kroupa_mass_distribution(1, in_primary_mass_max)[0]
            else: # flat distribution
               in_mass_ratio = flat_distr(in_mass_ratio_min, in_mass_ratio_max)
               self.in_secondary_mass = in_mass_ratio * self.in_primary_mass        


        if out_mass_ratio_max == out_mass_ratio_min:
            out_mass_ratio = out_mass_ratio_min
            self.out_mass = out_mass_ratio * (self.in_primary_mass + self.in_secondary_mass)
        else: 
            if out_mass_ratio_distr == 1:# Kroupa 2001 
                self.out_mass = new_kroupa_mass_distribution(1, in_primary_mass_max)[0]
#                while self.out_mass < in_primary_mass_min:
#                        self.out_mass = new_kroupa_mass_distribution(1, in_primary_mass_max)[0]
            else: # flat distribution
               out_mass_ratio = flat_distr(out_mass_ratio_min, out_mass_ratio_max)
               self.out_mass = out_mass_ratio * (self.in_primary_mass + self.in_secondary_mass)


    def generate_semi(self, 
                        in_semi_max, in_semi_min, 
                        out_semi_max, out_semi_min,
                        in_semi_distr,  out_semi_distr):
                        
        if in_semi_max == in_semi_min:
            self.in_semi = in_semi_min
        else:
            if in_semi_distr == 1: #Constant 
                if REPORT_USER_WARNINGS:
                    print 'TPS::generate_semi: unambiguous choise of constant semi-major axis'
                    print '--A_min option to set the value of the semi-major axis in the inner binary'                
                self.in_semi = in_semi_min
            elif in_semi_distr == 2: #Tokovinin
                self.in_semi = 0.|units.RSun
                while (self.in_semi < in_semi_min or self.in_semi > in_semi_max):
                    logP = np.random.normal(5, 2.3, 1)
                    P = (10**logP)|units.day
                    self.in_semi = ((P/2/np.pi)**2 * constants.G* (self.in_primary_mass + self.in_secondary_mass))**(1./3.)  
                    if logP < -0.3 or logP > 10:#truncation of Gaussian wings
                        self.in_semi = 0.|units.RSun
            else: # log flat distribution
                 maximal_semi = min(in_semi_max, out_semi_max)
                 self.in_semi = log_flat_distr(in_semi_min, maximal_semi)
                        
        
                        
        if out_semi_max == out_semi_min:
            self.out_semi = out_semi_min
        else:
            if out_semi_distr == 1: #Constant 
                if REPORT_USER_WARNINGS:
                    print 'TPS::generate_semi: unambiguous choise of constant semi-major axis'
                    print '--a_min option to set the value of the semi-major axis in the outer binary'                
                self.out_semi = out_semi_min
            elif out_semi_distr == 2: #Tokovinin
                self.out_semi = 0.|units.RSun
                while (self.out_semi < out_semi_min or self.out_semi > in_semi_max):
                    logP_out = np.random.normal(5, 2.3, 1)
                    P_out = (10**logP_out)|units.day
                    self.out_semi = ((P_out/2/np.pi)**2 * constants.G* (self.in_primary_mass + self.in_secondary_mass + self.out_mass))**(1./3.)                    
                    if logP_out < -0.3 or logP_out > 10:#truncation of Gaussian wings
                        self.out_semi = 0.|units.RSun
                    if logP_out < 3: # no bifurcation
                        self.out_semi = 0.|units.RSun

                #note:overwrites the inner orbital separation
                self.in_semi = 0.|units.RSun
                while (self.in_semi < in_semi_min or self.in_semi > in_semi_max):
                    logP_in = np.random.normal(5, 2.3, 1)
                    P_in = (10**logP_in)|units.day
                    self.in_semi = ((P_in/2/np.pi)**2 * constants.G* (self.in_primary_mass + self.in_secondary_mass))**(1./3.)                    
                    if logP_in < -0.3 or logP_in > 10:#truncation of Gaussian wings
                        self.in_semi = 0.|units.RSun

                    dlogP = logP_out - logP_in #unstable
                    if dlogP < 0.7:
                        self.in_semi = 0.|units.RSun   
                    elif dlogP < 1.7:
                        x = np.random.uniform(0, 1, 1)
                        if dlogP - 0.7 > x:
                            self.in_semi = 0.|units.RSun   
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
                if REPORT_USER_WARNINGS:
                    print 'TPS::generate_ecc: unambiguous choise of constant eccentricity'
                    print '--E_min option to set the value of the eccentricity in the inner binary'                
                self.in_ecc = in_ecc_min
            else: #Thermal distribution
                 self.in_ecc = np.sqrt(np.random.uniform(in_ecc_min, in_ecc_max))
                 
        if out_ecc_max == out_ecc_min:
            self.out_ecc = out_ecc_min
        else:
            if out_ecc_distr == 1: #Constant 
                if REPORT_USER_WARNINGS:
                    print 'TPS::generate_ecc: unambiguous choise of constant eccentricity'
                    print '--e_min option to set the value of eccentricity in the outer binary'                
                self.out_ecc = out_ecc_min
            else: #Thermal distribution
                 self.out_ecc = np.sqrt(np.random.uniform(out_ecc_min, out_ecc_max))


    def generate_incl(self, incl_max, incl_min, incl_distr):
        if incl_max == incl_min:
            self.incl = incl_min           
        else:
            if incl_distr == 1: #Constant 
                if REPORT_USER_WARNINGS:
                    print 'TPS::generate_incl: unambiguous choise of constant relative inclination'
                    print '--i_min option to set the value of the relative inclination in the inner triple'                
                self.incl = incl_min
            else: #Circular uniform distribution 
                 self.incl = np.arccos(np.random.uniform(np.cos(incl_min), np.cos(incl_max)))
                 
    def generate_aop(self,
                        in_aop_max, in_aop_min, 
                        out_aop_max, out_aop_min,
                        in_aop_distr, out_aop_distr):

        if in_aop_max == in_aop_min:
            self.in_aop = in_aop_min
        else:
            if in_aop_distr == 1: #Constant 
                if REPORT_USER_WARNINGS:
                    print 'TPS::generate_aop: unambiguous choise of constant argument of pericenter'
                    print '--G_min option to set the value of the argument of pericenter of the inner binary'                
                self.in_aop = in_aop_min
            else: #Uniform distribution 
                 self.in_aop = np.random.uniform(in_aop_min, in_aop_max)
                     
                 
        if out_aop_max == out_aop_min:
            self.out_aop = out_aop_min
        else:
            if out_aop_distr == 1: #Constant 
                if REPORT_USER_WARNINGS:
                    print 'TPS::generate_aop: unambiguous choise of constant argument of pericenter'
                    print '--g_min option to set the value of the argument of pericenter of the outer binary'                
                self.out_aop = out_aop_min
            else: #Uniform distribution 
                 self.out_aop = np.random.uniform(out_aop_min, out_aop_max)
                

    def generate_loan(self,
                        in_loan_max, in_loan_min, in_loan_distr):                                

        if in_loan_max == in_loan_min:
            self.in_loan = in_loan_min
        else:
            if in_loan_distr == 0: #Circular uniform distribution
                self.in_loan = np.arccos(np.random.uniform(np.cos(in_loan_min), np.cos(in_loan_max)))
            else: #Constant
                if REPORT_USER_WARNINGS:
                    print 'TPS::generate_loan: unambiguous choise of constant longitude of ascending nodes'
                    print '--O_min option to set the value of the argument of pericenter of the inner binary'                
                self.in_loan = in_loan_min
                    

#-------

# Eggleton 2009, 399, 1471
    def generate_mass_and_semi_eggleton(self, in_primary_mass_max, in_primary_mass_min, in_semi_max, in_semi_min, 
                        out_semi_max, out_semi_min):
#        U0_mass = [0., .01, .09, .32, 1., 3.2, 11, 32, np.inf]|units.MSun
#        U0_l0 = [0.40, 0.40, 0.40, 0.40, 0.50, 0.75, 0.88, 0.94, 0.96]#        U0_l1 = [0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.20, 0.60, 0.80]#        U0_l2 = [0.00, 0.00, 0.00, 0.00, 0.00, 0.20, 0.33, 0.82, 0.90]#        U0_l3 = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
#              
#        Mt = eggleton_mass_distr(0.1|units.MSun, 100|units.MSun)
#        U = np.random.uniform(0, 1)
#        f_l0 = interp1d (U0_mass, U0_l0)
#        U0 = f_l0(Mt)     

        U0_mass = [0., .01, .09, .32, 1., 3.2, 11, 32, np.inf]#solar mass
        U0_l0 = [0.40, 0.40, 0.40, 0.40, 0.50, 0.75, 0.88, 0.94, 0.96]        U0_l1 = [0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.20, 0.60, 0.80]        U0_l2 = [0.00, 0.00, 0.00, 0.00, 0.00, 0.20, 0.33, 0.82, 0.90]        U0_l3 = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
              
        Mt = eggleton_mass_distr(0.1|units.MSun, 100|units.MSun)
        U = np.random.uniform(0, 1)
        f_l0 = interp1d (U0_mass, U0_l0)
        U0 = f_l0(Mt.value_in(units.MSun)) #messy, but otherwise cluster crashes    

        while U >= U0:
            U = np.random.uniform(0, 1)
            f_l0 = interp1d (U0_mass, U0_l0)
            U0 = f_l0(Mt.value_in(units.MSun))     
        
        
        
        if U < U0:
            V = np.random.uniform(0, 1)
            P0 = 1.e5 * V**2 / (1-V)**2.5 |units.day
            while P0 > 1e10|units.day:             
                V = np.random.uniform(0, 1)
                P0 = 1.e5 * V**2 / (1-V)**2.5|units.day

            x_p = np.random.uniform(0, 1)
            if P0 > 25|units.day or x_p > 0.25: 
                Q0 = ( (U0-U)/U0 )**0.8
            else:
                Q0 = 0.9+0.09*(U0-U)/U0 
            if Q0 < 0.01:
                Q0 = 0.01                                
                
            M1 = Mt / (1+Q0)
            M2 = M1 * Q0
            f_l1 = interp1d(U0_mass, U0_l1)
            
            U1 = np.random.uniform(0, 1)
            U1_0 = f_l1(M1.value_in(units.MSun))     
            U2 = np.random.uniform(0, 1)
            U2_0 = f_l1(M2.value_in(units.MSun))
            
            #M1 bifurcutas and M2 not
            if U1< U1_0 and U2>=U2_0:
                M_bin = M1  
                U_bin = U1 
                U0_bin = U1_0
                M_comp = M2
            #M2 bifurcutas and M1 not
            elif U1>= U1_0 and U2<U2_0:
                M_bin = M2
                U_bin = U2
                U0_bin = U2_0
                M_comp = M1
            else: #two bifurcations -> higher order multiplicity
#                print U1, U1_0, U2, U2_0
#                exit(1)
                 return False
                
            V_bin = np.random.uniform(0, 1)
            P_bin = 0.2 * P0 * 10**(-5*V_bin)

            x_pb = np.random.uniform(0, 1)
            if P_bin > 25|units.day or x_pb > 0.25: 
                Q_bin = ( (U0_bin-U_bin)/U0_bin )**0.8
            else:
                Q_bin = 0.9+0.09*(U0_bin-U_bin)/U0_bin 
            if Q_bin < 0.01:
                Q_bin = 0.01                                

            M1_bin = M_bin / (1+Q_bin)
            M2_bin = M1_bin * Q_bin

            self.in_primary_mass = M1_bin
            self.in_secondary_mass = M2_bin
            self.out_mass = M_comp
            
            self.in_semi = ((P_bin/2*np.pi)**2 * M_bin*constants.G ) ** (1./3.)
            self.out_semi = ((P0/2*np.pi)**2 * Mt*constants.G ) ** (1./3.)
      
            if self.in_primary_mass < in_primary_mass_min or self.in_primary_mass > in_primary_mass_max:
                return False 
            if self.in_semi < in_semi_min or self.in_semi > in_semi_max:
                return False
            if self.out_semi < out_semi_min or self.out_semi > out_semi_max:
                return False
                                                         
            return True
        else:
            print 'not possible'
            exit(1)

                 

#-------
        
#-------
    def print_triple(self):
        print '\nTriple - ' 
        print 'm =', self.in_primary_mass, self.in_secondary_mass, self.out_mass
        print 'a =', self.in_semi, self.out_semi
        print 'e =', self.in_ecc, self.out_ecc
        print 'i =', self.incl
        print 'g =', self.in_aop, self.out_aop
        print 'o =', self.in_loan, self.in_loan -np.pi
#-------

#-------
def evolve_model(in_primary_mass_max, in_primary_mass_min, 
                        in_mass_ratio_max, in_mass_ratio_min, 
                        out_mass_ratio_max, out_mass_ratio_min, 
                        in_semi_max, in_semi_min, out_semi_max, out_semi_min, 
                        in_ecc_max, in_ecc_min, out_ecc_max, out_ecc_min, 
                        incl_max, incl_min,
                        in_aop_max, in_aop_min, out_aop_max, out_aop_min,
                        in_loan_max, in_loan_min, 
                        in_primary_mass_distr, in_mass_ratio_distr, out_mass_ratio_distr,
                        in_semi_distr,  out_semi_distr, in_ecc_distr, out_ecc_distr, incl_distr,
                        in_aop_distr, out_aop_distr, in_loan_distr,                                                                      
                        metallicity, tend, number, initial_number, seed,
                        stop_at_merger, stop_at_disintegrated, stop_at_triple_mass_transfer,
                        stop_at_collision, stop_at_dynamical_instability, stop_at_mass_transfer,
                        file_name, file_type):


    
    i_n = 0
    nr_ids = 0 #number of systems that is dynamically unstable at initialisation
    nr_imt = 0 #number of systems that has mass transfer at initialisation
    while i_n < number:
        triple_system = Generate_initial_triple(in_primary_mass_max, in_primary_mass_min, 
                    in_mass_ratio_max, in_mass_ratio_min, 
                    out_mass_ratio_max, out_mass_ratio_min, 
                    in_semi_max, in_semi_min, out_semi_max, out_semi_min, 
                    in_ecc_max, in_ecc_min, out_ecc_max, out_ecc_min, 
                    incl_max, incl_min,
                    in_aop_max, in_aop_min, out_aop_max, out_aop_min,
                    in_loan_max, in_loan_min, 
                    in_primary_mass_distr, in_mass_ratio_distr, out_mass_ratio_distr,
                    in_semi_distr,  out_semi_distr, in_ecc_distr, out_ecc_distr, incl_distr,
                    in_aop_distr, out_aop_distr, in_loan_distr)
                
        if REPORT:
           triple_system.print_triple()

        if (min_mass > triple_system.in_primary_mass) or (min_mass > triple_system.in_secondary_mass) or (min_mass > triple_system.out_mass):
                if REPORT:
                    print 'non-star included: ', triple_system.in_primary_mass, triple_system.in_secondary_mass, triple_system.out_mass
                continue

        number_of_system = initial_number + i_n
        if REPORT:
            print 'number of system = ', number_of_system
        tr = triple.main(inner_primary_mass = triple_system.in_primary_mass, 
                    inner_secondary_mass = triple_system.in_secondary_mass, 
                    outer_mass = triple_system.out_mass, 
                    inner_semimajor_axis = triple_system.in_semi, 
                    outer_semimajor_axis = triple_system.out_semi, 
                    inner_eccentricity = triple_system.in_ecc, 
                    outer_eccentricity = triple_system.out_ecc, 
                    relative_inclination = triple_system.incl, 
                    inner_argument_of_pericenter = triple_system.in_aop, 
                    outer_argument_of_pericenter = triple_system.out_aop, 
                    inner_longitude_of_ascending_node = triple_system.in_loan, 
                    metallicity = metallicity, tend = tend, number = number_of_system, 
                    stop_at_merger = stop_at_merger, stop_at_disintegrated = stop_at_disintegrated,
                    stop_at_triple_mass_transfer = stop_at_triple_mass_transfer, stop_at_collision = stop_at_collision, 
                    stop_at_dynamical_instability = stop_at_dynamical_instability, stop_at_mass_transfer = stop_at_mass_transfer, 
                    stop_at_init_mass_transfer = stop_at_init_mass_transfer,
                    file_name = file_name, file_type = file_type)                        

        if tr.triple.dynamical_instability_at_initialisation == True:
            nr_ids +=1
        elif tr.triple.mass_transfer_at_initialisation == True:
            nr_imt +=1
        else:
            i_n += 1            

    if REPORT:
      print number, i_n, nr_ids, nr_imt                               

def test_initial_parameters(in_primary_mass_max, in_primary_mass_min, 
                        in_mass_ratio_max, in_mass_ratio_min, 
                        out_mass_ratio_max, out_mass_ratio_min, 
                        in_semi_max, in_semi_min, out_semi_max, out_semi_min, 
                        in_ecc_max, in_ecc_min,  out_ecc_max, out_ecc_min, 
                        incl_max, incl_min,
                        in_aop_max, in_aop_min, out_aop_max, out_aop_min,
                        in_loan_max, in_loan_min, 
                        in_primary_mass_distr, in_mass_ratio_distr, out_mass_ratio_distr,
                        in_semi_distr,  out_semi_distr, in_ecc_distr, out_ecc_distr, incl_distr,
                        in_aop_distr, out_aop_distr, in_loan_distr,                    
                        metallicity, tend, number, initial_number, seed,
                        stop_at_merger, stop_at_disintegrated, stop_at_triple_mass_transfer,
                        stop_at_collision, stop_at_dynamical_instability, stop_at_mass_transfer,
                        file_name, file_type):

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
        print 'error: relative inclination not in allowed range [0, pi]'
        exit(1)

    if (incl_max < incl_min):
        print 'error: maximum relative inclination smaller than minimum relative inclination'
        exit(1)



    if (in_aop_min < -np.pi) or (in_aop_max > np.pi):
        print 'error: inner argument of pericenter not in allowed range [-pi,pi]'
        exit(1)

    if (out_aop_min < -np.pi) or (out_aop_max > np.pi):
        print 'error: outer argument of pericenter not in allowed range [-pi,pi]'
        exit(1)

    if (in_aop_max < in_aop_min):
        print 'error: maximum inner argument of pericenter smaller than minimum argument of pericenter'
        exit(1)

    if (out_aop_max < out_aop_min):
        print 'error: maximum outer argument of pericenter smaller than minimum argument of pericenter'
        exit(1)


    if (in_loan_min < -1*np.pi) or (in_loan_max > np.pi):
        print 'error: inner longitude of ascending node not in allowed range [-pi,pi]'
        exit(1)

    if (in_loan_max < in_loan_min):
        print 'error: maximum inner longitude of ascending node smaller than minimum argument of pericenter'
        exit(1)

    if (number < 1):
        print 'Requested number of systems < 1'
        exit(1)

    if (initial_number < 0):
        print 'Initial number of system < 0'
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
                      
    parser.add_option("--Q_max", dest="in_mass_ratio_max", type="float", default = 1.0,
                      help="maximum of inner mass ratio [%default]")
    parser.add_option("--Q_min", dest="in_mass_ratio_min", type="float", default = 0.,
                      help="minimum of inner mass ratio [%default]")
    parser.add_option("--Q_distr", dest="in_mass_ratio_distr", type="int", default = 0,
                      help="inner mass ratio distribution [Flat]")

    parser.add_option("--q_max", dest="out_mass_ratio_max", type="float", default = 1.0,
                      help="maximum of outer mass ratio [%default]")
    parser.add_option("--q_min", dest="out_mass_ratio_min", type="float", default = 0.,
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
                      help="inner semimajor axis distribution [logFlat]")


    parser.add_option("--a_min", unit=units.RSun,
                      dest="out_semi_min", type="float", 
                      default = 5|units.RSun,
                      help="minimum of outer semi major axis [%default]")
    parser.add_option("--a_max", unit=units.RSun,
                      dest="out_semi_max", type="float", 
                      default = 5e6|units.RSun,
                      help="maximum of outer semi major axis [%default]")
    parser.add_option("--a_distr", dest="out_semi_distr", type="int", default = 0,
                      help="outer semimajor axis distribution [logFlat]")

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
                      help="minimum of relative inclination [rad] [%default]")
    parser.add_option("--i_max, --I_max",
                      dest="incl_max", type="float", default = np.pi,
                      help="maximum of relative inclination [rad] [%default]")
    parser.add_option("--i_distr, --I_distr", dest="incl_distr", type="int", default = 0,
                      help="relative inclination distribution [Circular uniform]")

                      
    parser.add_option("--G_min",
                      dest="in_aop_min", type="float", default = -np.pi,
                      help="minimum of inner argument of pericenter [rad] [%default]")
    parser.add_option("--G_max",
                      dest="in_aop_max", type="float", default = np.pi,
                      help="maximum of inner argument of pericenter [rad] [%default]")
    parser.add_option("--G_distr", dest="in_aop_distr", type="int", default = 0,
                      help="inner argument of pericenter distribution [Uniform]")

    parser.add_option("--g_min",
                      dest="out_aop_min", type="float", default = -np.pi,
                      help="minimum of outer argument of pericenter [rad] [%default]")
    parser.add_option("--g_max",
                      dest="out_aop_max", type="float", default = np.pi,
                      help="maximum of outer argument of pericenter [rad] [%default]")
    parser.add_option("--g_distr", dest="out_aop_distr", type="int", default = 0,
                      help="outer argument of pericenter distribution [Uniform]")
                      
    parser.add_option("--O_min",
                      dest="in_loan_min", type="float", default = -np.pi,
                      help="minimum of inner longitude of ascending node [rad] [%default]")
    parser.add_option("--O_max",
                      dest="in_loan_max", type="float", default = np.pi,
                      help="maximum of inner longitude of ascending node [rad] [%default]")
    parser.add_option("--O_distr", dest="in_loan_distr", type="int", default = 1,
                      help="inner longitude of ascending node distribution [Constant]")
                     
    parser.add_option("-t", "-T", unit=units.Myr, 
                      dest="tend", type="float", default = 13500|units.Myr,
                      help="end time [%default] %unit")
    parser.add_option("-z", unit=units.none, 
                      dest="metallicity", type="float", default = 0.02|units.none,
                      help="metallicity [%default] %unit")
    parser.add_option("-n", dest="number", type="int", default = 10,
                      help="number of systems [%default]")
    parser.add_option("-N", dest="initial_number", type="int", default = 0,
                      help="number of initial system [%default]")
    parser.add_option("-s", unit=units.none, 
                      dest="seed", type="float", default = 0.|units.none,
                      help="seed [%default] %unit")
#    int actual_seed = srandinter(input_seed);
    parser.add_option("--stop_at_merger", dest="stop_at_merger", action="store_false", default = True, 
                      help="stop at merger [%default] %unit")
    parser.add_option("--stop_at_disintegrated", dest="stop_at_disintegrated", action="store_false", default = True,
                      help="stop at disintegrated [%default] %unit")
    parser.add_option("--stop_at_triple_mass_transfer", dest="stop_at_triple_mass_transfer", action="store_false", default = True,
                      help="stop at triple mass transfer [%default] %unit")
    parser.add_option("--stop_at_collision", dest="stop_at_collision", action="store_false",default = True,
                      help="stop at collision [%default] %unit")
    parser.add_option("--stop_at_dynamical_instability", dest="stop_at_dynamical_instability", action="store_false", default = True,
                      help="stop at dynamical instability [%default] %unit")
    parser.add_option("--stop_at_mass_transfer", dest="stop_at_mass_transfer", action="store_true", default = False,
                      help="stop at mass transfer [%default] %unit")
    parser.add_option("-f", dest="file_name", type ="string", default = "triple_nokozaidt.hdf",#"triple.txt"
                      help="file name[%default]")
    parser.add_option("-F", dest="file_type", type ="string", default = "hdf5",#"txt"
                      help="file type[%default]")



    options, args = parser.parse_args()
    return options.__dict__



if __name__ == '__main__':
    options = parse_arguments()
    set_printing_strategy("custom", 
                          preferred_units = [units.MSun, units.RSun, units.Myr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")

    test_initial_parameters(**options)
    evolve_model(**options)



