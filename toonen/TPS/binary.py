from amuse.lab import *
import numpy as np
from math import sqrt

REPORT_BINARY_EVOLUTION = False
REPORT_FUNCTION_NAMES = False

#constants
which_common_envelope = 0
common_envelope_efficiency = 1.0
envelope_structure_parameter = 0.5
common_envelope_efficiency_gamma = 1.75

stellar_types_giants = [2,3,4,5,6,8,9]
stellar_types_remnants = [10,11,12,13,14]


def print_star(star):
    print star.age, star.stellar_type, star.mass, star.radius, star.core_mass, star.spin_angular_frequency,

def print_binary(bs):
    print binary.semi_major_axis, binary.eccentricity, binary.argument_of_pericenter, binary.longitude_of_ascending_node
    print '\t'
    print_star(bs.child1)
    print '\t'
    print_star(bs.child2)
    print '\n'
    
    
def print_triple(tr):
    print tr.mutual_inclination
    print_binary(tr.inner_binary)
    print_binary(tr.outer_binary)
    
def roche_radius_dimensionless(M, m) :
    # Assure that the q is calculated in identical units.
    unit = M.unit
    # and that q itself has no unit
    q = M.value_in(unit)/m.value_in(unit)
    q13 =  q**(1./3.)
    q23 =  q13**2
    return  0.49*q23/(0.6*q23 + np.log(1 + q13))

def roche_radius(bin, primary):
    if bin.is_binary and primary.is_star:
        return bin.semimajor_axis * roche_radius_dimensionless(primary.mass, bin.mass-primary.mass)

    print 'Error: Roche radius can only be determined in a binary'
    exit(1)


def common_envelope_efficiency(donor, accretor):
    return common_envelope_efficiency

def envelope_structure_parameter(donor):
    return envelope_structure_parameter
    
def common_envelope_efficiency_gamma(donor, accretor):
    return common_envelope_efficiency_gamma
    

# ang.mom balance: \Delta J = \gamma * J * \Delta M / M
# See Eq. 5 of Nelemans VYPZ 2000, 360, 1011 A&A
def common_envelope_angular_momentum_balance(bs, donor, accretor, triple):
    if REPORT_FUNCTION_NAMES:
        print 'Common envelope angular momentum balance'

    if REPORT_BINARY_EVOLUTION:
        if bs.eccentricity > 0.05:
            print 'gamma common envelope in eccentric binary'
            exit(-1)
        print 'Before common envelope angular momentum balance' 
        print_binary(bs) 

    gamma = common_envelope_efficiency_gamma(donor, accretor)
    J_init = sqrt(bs.semi_major_axis) * (donor.mass * accretor.mass) / sqrt(donor.mass + accretor.mass) * sqrt(1-bs.eccentricity**2)
    J_f_over_sqrt_a_new = (donor.core_mass * accretor.mass) / sqrt(donor.core_mass + accretor.mass)
    J_lost = gamma * donor.envelope_mass * J_i/(donor.mass + accretor.mass)
    sqrt_a_new = max(0, (J_i -J_lost)/J_f_over_sqrt_a_new)
    a_new = pow(sqrt_a_new, 2)

    Rl_donor_new = roche_radius_dimensionless(donor.core_mass, accretor.mass)*a_new
    Rl_accretor_new = roche_radius_dimensionless(accretor.mass, donor.core_mass)*a_new    
    
    if (donor.core_radius > Rl_donor_new) or (accretor.radius > Rl_accretor_new):
        print 'Merger in inner binary through common envelope phase (energy balance)'
        exit(0)
    else:
        bs.semi_major_axis = a_new
        bs.eccentricity = zero
#        bs.argument_of_pericenter = 
#        bs.inner_longitude_of_ascending_node =  
#        donor.spin_angular_frequency = 
#        accretor.spin_angular_frequency = 
        donor_in_se_code = donor.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
        #donor_in_se_code.change_mass(-1*donor.envelope_mass)    reduce_mass not subtrac mass     
        
        triple.channel_from_se.copy()
        if REPORT_BINARY_EVOLUTION:
            print 'After common envelope angular momentum balance' 
            print_binary(bs) 
        donor.is_donor = False
        


#Following Webbink 1984
def common_envelope_energy_balance(bs, donor, accretor, triple):
    if REPORT_FUNCTION_NAMES:
        print 'Common envelope energy balance'

    if REPORT_BINARY_EVOLUTION:
        print 'Before common envelope angular momentum balance' 
        print_binary(bs) 

    alpha = common_envelope_efficiency(donor, accretor) 
    lambda_donor = envelope_structure_parameter(donor)

    Rl_donor = roche_radius(bs, donor)
    donor_radius = min(donor.radius, Rl_donor)

    orb_energy_new = donor.mass * donor.envelope_mass / (alpha * lambda_donor * donor_radius) + donor.mass * accretor.mass/2/bs.semi_major_axis
    a_new = donor.core_mass * accretor.mass / 2 / orb_energy_new
#    a_new = bs.semi_major_axis * (donor.core_mass/donor.mass) / (1. + (2.*donor.envelope_mass*bs.semi_major_axis/(alpha_lambda*donor_radius*accretor.mass)))
    
    Rl_donor_new = roche_radius_dimensionless(donor.core_mass, accretor.mass)*a_new
    Rl_accretor_new = roche_radius_dimensionless(accretor.mass, donor.core_mass)*a_new    
    
    if (donor.core_radius > Rl_donor_new) or (accretor.radius > Rl_accretor_new):
        print 'Merger in inner binary through common envelope phase (energy balance)'
        exit(0)
    else:
        bs.semi_major_axis = a_new
        bs.eccentricity = zero
#        bs.argument_of_pericenter = 
#        bs.inner_longitude_of_ascending_node =  
#        donor.spin_angular_frequency = 
#        accretor.spin_angular_frequency = 
        donor_in_se_code = donor.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
        #donor_in_se_code.change_mass(-1*donor.envelope_mass)    reduce_mass not subtrac mass     
        
        triple.channel_from_se.copy()
        if REPORT_BINARY_EVOLUTION:
            print 'After common envelope energy balance' 
            print_binary(bs) 

        donor.is_donor = False


# See appendix of Nelemans YPZV 2001, 365, 491 A&A
def double_common_envelope_energy_balance(bs, donor, accretor, triple):
    if REPORT_FUNCTION_NAMES:
        print 'Double common envelope energy balance'

    if REPORT_BINARY_EVOLUTION:
        print 'Before common envelope angular momentum balance' 
        print_binary(bs) 

    alpha = common_envelope_efficiency(donor, accretor)
    lambda_donor = envelope_structure_parameter(donor) 
    lambda_accretor = envelope_structure_parameter(accretor)

    Rl_donor = roche_radius(bs, donor)
    donor_radius = min(donor.radius, Rl_donor)
    accretor_radius = accretor.radius
    
    orb_energy_new = donor.mass * donor.envelope_mass / (alpha * lambda_donor * donor_radius) + accretor.mass * accretor.envelope_mass / (alpha * lambda_accretor * accretor_radius) + donor.mass * accretor.mass/2/bs.semi_major_axis
    a_new = donor.core_mass * accretor.core_mass / 2 / orb_energy_new

    Rl_donor_new = roche_radius_dimensionless(donor.core_mass, accretor.core_mass)*a_new
    Rl_accretor_new = roche_radius_dimensionless(accretor.core_mass, donor.core_mass)*a_new    
    
    if (donor.core_radius > Rl_donor_new) or (accretor.core_radius > Rl_accretor_new):
        print 'Merger in inner binary through common envelope phase (double common envelope)'
        exit(0)
    else:
        bs.semi_major_axis = a_new
        bs.eccentricity = zero
#        bs.argument_of_pericenter = 
#        bs.inner_longitude_of_ascending_node =  
#        donor.spin_angular_frequency = 
#        accretor.spin_angular_frequency = 
        donor_in_se_code = donor.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
        #donor_in_se_code.change_mass(-1*donor.envelope_mass)    reduce_mass not subtrac mass     
        accretor_in_se_code = accretor.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
        #accretor_in_se_code.change_mass(-1*accretor.envelope_mass)    reduce_mass not subtrac mass     

        triple.channel_from_se.copy()
        if REPORT_BINARY_EVOLUTION:
            print 'After double common envelope energy balance' 
            print_binary(bs) 

        donor.is_donor = False



def common_envelope_phase(bs, donor, accretor, triple):
    if REPORT_FUNCTION_NAMES:
        print 'Common envelope phase'

    if donor.stellar_type not in stellar_types_giants and accretor.stellar_type not in stellar_types_giants:
#        possible options: MS+MS, MS+remnant, remnant+remnant,
#                          HeMS+HeMS, HeMS+MS, HeMS+remnant
        print donor.stellar_type, accretor_stellar_type
        print 'Merger in inner binary through common envelope phase (stellar types)'
        exit(0)
        
    
    if which_common_envelope == 0:
        if donor.stellar_type in stellar_types_giants and accretor.stellar_type in stellar_types_giants:
            double_common_envelope(bs, donor, accretor, triple)
        else:
            common_envelope_energy_balance(bs, donor, accretor, triple)
    elif which_common_envelope == 1:
        if donor.stellar_type in stellar_types_giants and accretor.stellar_type in stellar_types_giants:
            double_common_envelope(bs, donor, accretor, triple)
        else:
            common_envelope_angular_momentum_balance(bs, donor, accretor, triple)
    elif which_common_envelope == 2:
        if donor.stellar_type in stellar_types_giants and accretor.stellar_type in stellar_types_giants:
            #giant+giant
            double_common_envelope(bs, donor, accretor, triple)
        elif donor.stellar_type in stellar_types_remnants and accretor.stellar_type in stellar_types_remnants:
            #giant+remnant
            common_envelope_energy_balance(bs, donor, accretor, triple)
        else:
            #giant+normal(non-giant, non-remnant)
            common_envelope_angular_momentum_balance(bs, donor, accretor, triple)
    

def stable_mass_transfer(bs, donor, accretor, triple):
    # orbital evolution is being taken into account in secular_code        
    if REPORT_FUNCTION_NAMES:
        print 'Stable mass transfer'

    Md = donor.mass
    Ma = accretor.mass
    print Md, Ma, 
    print donor.previous_mass, accretor.previous_mass

    #set mass transfer rate
    
    dt = tr.timestep
    dm = bs.mass_transfer_rate * dt
    print 'check sign of dm:', dm, 'presumably dm >0'
    donor_in_se_code = donor.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
    donor_in_se_code.change_mass(-1*dm, dt)
    
    # there is an implicit assumption in change_mass that the accreted mass is of solar composition (hydrogen)
    accretor_in_se_code = accretor.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
#    accretor_in_se_code.change_mass(dm, dt)
    # for now, only conservative mass transfer   
    accretor_in_se_code.change_mass(dm, -1*dt)

    triple.channel_from_se.copy()

    Md_new = donor.mass
    Ma_new = accretor.mass
    print Md_new, Ma_new, Md-Md_new, Ma-Ma_new
    print 'check if this is indeed conservative'
    accretion_efficiency = (Ma_new-Ma)/(Md-Md_new)
    print accretion_efficiency
#    bs.accretion_efficiency_mass_transfer = accretion_efficiency



def orbital_angular_momentum(bs):

    M = bs.child1.mass
    m = bs.child2.mass
    a = bs.semimajor_axis
    e = bs.eccentricity
    J = M*m * np.sqrt(constants.G*a*(1-e**2)/(M+m))

    if REPORT_BINARY_EVOLUTION:
        print 'Jorb:', M, m, a, e, J

    return J

def stellar_angular_momentum(ss):

    moment_of_inertia = ss.gyration_radius**2 * ss.mass * ss.radius**2
    Jstar = moment_of_inertia * ss.spin_angular_frequency
    return Jstar


def orbital_period(bs):

    Porb = 2*np.pi * np.sqrt(bs.semimajor_axis**3/constants.G / (bs.child1.mass + bs.child2.mass))
    return Porb


def mass_transfer_stability_binary(bs, donor, accretor, tr):

    if REPORT_FUNCTION_NAMES:
        print 'Mass transfer stability'

    Js = stellar_angular_momentum(donor)#accretor ook? check this
    Jb = orbital_angular_momentum(bs)
    
    if REPORT_BINARY_EVOLUTION:
        print "Darwin Riemann instability?:", Js, Jb, Jb/3.
    
        
    #which outer_bs.child is bs..?
    outer_bs = tr.outer_binary # niet mooi
    if outer_bs.child1.is_star and outer_bs.child1.is_donor:
        print 'RLOF in inner and outer binary'
        exit(0)
    if outer_bs.child2.is_star and outer_bs.child2.is_donor:
        print 'RLOF in inner and outer binary'
        exit(0)
        
    
    if Js >= Jb/3. :
        if REPORT_BINARY_EVOLUTION:
            print "Darwin Riemann instability"
        common_envelope(bs, donor, accretor, tr)
        # heb ik zeker een triple?
        adjust_triple_after_ce_in_inner_binary(outer_bs, bs, tr)
        
    else :
        stable_mass_transfer(bs, donor, accretor, tr)
        #adjusting triple is done in detached


def adjust_triple_after_ce_in_inner_binary(bs, ce_binary, tr):
# Assumption:
# Unstable mass transfer (common-envelope phase) in the inner binary, affects the outer binary as a wind. 
# Instanteneous effect
   
    if REPORT_FUNCTION_NAMES:
        print 'Adjust triple after ce in inner_binary'


    M_com_after_ce = ce_binary.mass
    M_com_before_ce = ce_binary.previous_mass
    print M_com_after_ce, M_com_before_ce
    
    # accretion_efficiency
    M_accretor_before_ce = bs.child1.mass #niet zo mooi
    M_accretor_after_ce = bs.child1.mass #niet zo mooi
    
    adiabatic_expansion_due_to_mass_loss(bs.semi_major_axis, M_com_after_ce, M_com_before_ce, M_accretor_after_ce, M_accretor_before_ce)
#    bs.eccentricity = zero
#    bs.argument_of_pericenter = 
#    bs.inner_longitude_of_ascending_node =  
#    bs.child1.spin_angular_frequency = 

    Rl1, Rl2, Rl3 = triple.secular_code.give_roche_radii(triple.particles[0])
    if bs.child1.radius > Rl3: # niet zo mooi
        bs.child1.is_donor = True
    

def adiabatic_expansion_due_to_mass_loss(a_i, Md_f, Md_i, Ma_f, Ma_i):

    d_Md = Md_f - Md_i #negative mass loss rate
    d_Ma = Ma_f - Ma_i #positive mass accretion rate
    print "Stellar wind primary:", d_Md, d_Ma    

    Mt_f = Md_f + Ma_f
    Mt_i = Md_i + Ma_i

    if d_Md < 0 and d_Ma >= 0:
        eta = d_Ma / d_Md
        a_f = a_i * ((Md_f/Md_i)**eta * (Ma_f/Ma_i))**-2 * Mt_i/Mt_f
        return a_f

    return a_i
   

def contact_binary():
    if REPORT_FUNCTION_NAMES:
        print "Contact binary"

    #which outer_bs.child is bs..?
    if outer_bs.child1.is_star and outer_bs.child1.is_donor:
        print 'RLOF in inner and outer binary'
        exit(0)
    if outer_bs.child2.is_star and outer_bs.child2.is_donor:
        print 'RLOF in inner and outer binary'
        exit(0)


    print "Contact binary "
    exit(0)




## Calculates stellar wind velocoty.
## Steller wind velocity is 2.5 times stellar escape velocity
#def wind_velocity(star):
#    v_esc2 = constants.G * star.mass / star.radius
#    return 2.5*sqrt(v_esc2)
#}
#
#
## Bondi, H., and Hoyle, F., 1944, MNRAS 104, 273 (wind accretion.
## Livio, M., Warner, B., 1984, The Observatory 104, 152.
#def accretion_efficiency_from_stellar_wind(accretor, donor):
#velocity needs to be determined -> velocity average?
# why is BH dependent on ecc as 1/sqrt(1-e**2)

#    alpha_wind = 0.5
#    v_wind = wind_velocity(donor)
#    acc_radius = (constants.G*accretor.mass)**2/v_wind**4
#    
#    wind_acc = alpha_wind/sqrt(1-bs.eccentricity**2) / bs.semi_major_axis**2
#    v_factor = 1/((1+(velocity/v_wind)**2)**3./2.)
#    mass_fraction = acc_radius*wind_acc*v_factor
#
#    print 'mass_fraction:', mass_fraction
##    mass_fraction = min(0.9, mass_fraction)
#



def detached(bs, tr):
    # orbital evolution is being taken into account in secular_code        
    if REPORT_FUNCTION_NAMES:
        print 'Detached'
    
    # wind mass loss is done by se_code
    # wind accretion here:
    # update accretion efficiency of wind mass loss
    if bs.child1.is_star and bs.child2.is_star:
        bs.accretion_efficiency_wind_child1_to_child2 = 0.0
        bs.accretion_efficiency_wind_child2_to_child1 = 0.0

#        child1_in_se_code = bs.child1.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
#        child2_in_se_code = bs.child2.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
#
#        dt = tr.timestep
#        dm_child1_to_child2 = -1 * child1.wind_mass_loss_rate * bs.accretion_efficiency_wind_child1_to_child2 * dt
#        child2_in_se_code.change_mass(dm_child1_to_child2, dt)
#        dm_child12to_child1 = -1 * child2.wind_mass_loss_rate * bs.accretion_efficiency_wind_child2_to_child1 * dt
#        child1_in_se_code.change_mass(dm_child2_to_child1, dt)
# check if this indeed is accreted conservatively        


    elif bs.child1.is_star and bs.child2.is_binary and bs.child2.child1.is_star and bs.child2.child2.is_star:
        #Assumption: an inner binary is not effected by wind from an outer star
        bs.accretion_efficiency_wind_child1_to_child2 = 0.0

        bs.accretion_efficiency_wind_child2_to_child1 = 0.0
        
#        child1_in_se_code = bs.child1.as_set().get_intersecting_subset_in(self.se_code.particles)[0]
#        dt = tr.timestep
        
         #effect of wind from bs.child2.child1 onto bs.child1
#        mtr_w_in1_1 =  bs.child2.child1.wind_mass_loss_rate * (1-bs.child2.accretion_efficiency_wind_child1_to_child2)       
#        beta_w_in1_1 = 0.0
#        dm_in1_1 = -1 * mtr_w_in1_1 * beta_w_in1_1 * dt
#        
         #effect of wind from bs.child2.child2 onto bs.child1
#        mtr_w_in2_1 =  bs.child2.child2.wind_mass_loss_rate * (1-bs.child2.accretion_efficiency_wind_child2_to_child1)       
#        beta_w_in2_1 = 0.0
#        dm_in2_1 = -1 * mtr_w_in2_1 * beta_w_in2_1 * dt
#                    
#        dm = dm_in1_1 + dm_in2_1  
#        mtr = mtr_w_in1_1 + mtr_w_in2_1)


         #effect of mass transfer in the binary bs.child2 onto bs.child1
#        if bs.child2.child1.is_donor and bs.child2.child2.is_donor:
#            print 'contact binary in detached...'
#            exit(-1)
#        elif bs.child2.child1.is_donor or bs.child2.child2.is_donor:
#            #Assumption:
#            #Stable mass transfer in the inner binary, affects the outer binary as a wind.
#            mtr_rlof_in_1 = bs.child2.mass_transfer_rate * (1-bs.child2.accretion_efficiency_mass_transfer)
#            beta_rlof_in_1 = 0.0
#            dm_rlof_in_1 = -1 * mtr_rlof_in_1 * beta_rlof_in_1 * dt
#            dm += dm_rlof_in_1
#            mtr += mtr_rlof_in_1 

#        bs.accretion_efficiency_wind_child2_to_child1 = dm / ( mtr* -1 * dt)
            
#        child1_in_se_code.change_mass(dm, dt)
# check if this indeed is accreted conservatively        

    else:
        print 'detached: type of system unknown'
        print bs.child1.is_binary, bs.child1.is_star
        print bs.child2.is_binary, bs.child2.is_star
        exit(-1)                    
        
        
    #reset parameters after mass transfer
    tr.first_contact = True
    tr.particles[0].inner_binary.child1.mass_transfer_rate = 0.0 | units.MSun/units.yr
    tr.particles[0].inner_binary.child2.mass_transfer_rate = 0.0 | units.MSun/units.yr
    tr.particles[0].outer_binary.child1.mass_transfer_rate = 0.0 | units.MSun/units.yr




def triple_mass_transfer():
    # orbital evolution is being taken into account in secular_code        
    if REPORT_FUNCTION_NAMES:
        print 'Triple mass transfer'

    print "Mass transfer in outer binary"
    exit(0)




def resolve_binary_interaction(bs, tr):
   if REPORT_FUNCTION_NAMES:
        print 'Resolve binary interaction'
    
   if bs.is_binary and bs.child1.is_star:
        if REPORT_BINARY_EVOLUTION:
            Rl1 = roche_radius(bs, bs.child1)
            print "Check for RLOF:", bs.child1.mass, bs.child1.previous_mass
            print "Check for RLOF:", Rl1, bs.child1.radius
                
        if bs.child2.is_star:
            if REPORT_BINARY_EVOLUTION:
                Rl2 = roche_radius(bs, bs.child2)
                print "Check for RLOF:", bs.child2.mass, bs.child2.previous_mass
                print "Check for RLOF:", Rl2, bs.child2.radius

            if bs.child1.is_donor and bs.child2.is_donor:
                contact_binary()
            elif bs.child1.is_donor and not bs.child2.is_donor:
                mass_transfer_stability_binary(bs, bs.child1, bs.child2, tr)
            elif not bs.child1.is_donor and bs.child2.is_donor:
                mass_transfer_stability_binary(bs, bs.child2, bs.child1, tr)
            else:
                detached(bs, tr)
                                        
        elif bs.child2.is_binary:
            if REPORT_BINARY_EVOLUTION:
                print bs.mass, bs.child1.mass, bs.child2.mass
    

            if bs.child1.is_donor:
                if bs.child2.child1.is_donor or bs.child2.child2.is_donor:
                    print 'rlof in inner and outer binary'
                    exit(-1)
    
                triple_mass_transfer()
            else:
                detached(bs, tr)
                
        else:
            print 'resolve binary interaction: type of system unknown'
            print bs.is_binary, 
            print bs.child1.is_binary, bs.child1.is_star, 
            print bs.child2.is_binary, bs.child2.is_star
            exit(-1) 
                               
   else:
        print 'resolve binary interaction: type of system unknown'
        print bs.is_binary, 
        print bs.child1.is_binary, bs.child1.is_star, bs.child1.is_donor
        exit(-1)                    
                        



