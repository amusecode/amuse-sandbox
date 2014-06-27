from amuse.lab import *
import numpy as np
from math import sqrt

REPORT_BINARY_EVOLUTION = False
REPORT_FUNCTION_NAMES = False

#constants
which_common_envelope = 0
common_envelope_efficiency = 1.0
envelope_structure_parameter = 0.5

stellar_types_giants = [2,3,4,5,6,8,9]
stellar_types_remnants = [10,11,12,13,14]


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

def common_envelope_angular_momentum_balance(bs, donor, accretor, triple):
    if REPORT_FUNCTION_NAMES:
        print 'Common envelope angular momentum balance'

    #inner binary
    #remove mass
    #update star / se code
    #change inner orbit

def common_envelope_energy_balance(bs, donor, accretor, triple):
    if REPORT_FUNCTION_NAMES:
        print 'Common envelope energy balance'

    alpha_lambda = common_envelope_efficiency(donor, accretor) * envelope_structure_parameter(donor)
    Rl_donor = roche_radius(bs, donor)
    radius = min(donor.radius, Rl_donor)
    a_new = bs.semi_major_axis * (donor.core_radius/donor.mass) / (1. + (2.*donor.envelope_mass/(alpha_lambda*radius*accretor.mass)))
    
    Rl_donor_new = roche_radius_dimensionless(donor.mass, accretor.mass)*a_new
    Rl_accretor_new = roche_radius_dimensionless(accretor.mass, donor.mass)*a_new    
    
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



def double_common_envelope_energy_balance(bs, donor, accretor, triple):
    if REPORT_FUNCTION_NAMES:
        print 'Double common envelope energy balance'

    #inner binary
    #remove mass
    #update star / se code
    #change inner orbit

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


        
    
    #outer binary
    #adiabatic_expansion_due_to_mass_loss ->instantaneous effect

    donor.is_donor = False
    

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


def mass_transfer_stability(bs, donor, accretor, tr):

    if REPORT_FUNCTION_NAMES:
        print 'Mass transfer stability'

    Js = stellar_angular_momentum(donor)#accretor ook? check this
    Jb = orbital_angular_momentum(bs)
    
    if REPORT_BINARY_EVOLUTION:
        print "Darwin Riemann instability?:", Js, Jb, Jb/3.
    
    if Js >= Jb/3. :
        if REPORT_BINARY_EVOLUTION:
            print "Darwin Riemann instability"
        common_envelope(bs, donor, accretor, tr)
    else :
        stable_mass_transfer(bs, donor, accretor, tr)
#        accretor = companion_star(bs, donor)
#        rochelobe_overflow(bs, donor, accretor)




def semi_detached(bs, donor, accretor, tr):
    if REPORT_FUNCTION_NAMES:
        print 'Semi detached'

    print "Roche-lobe overflow by star with mass:", donor.mass
    mass_transfer_stability(bs, donor, accretor, tr)

    exit(0)



def contact_binary():
    if REPORT_FUNCTION_NAMES:
        print "Contact binary"

    print "Contact binary "
    exit(0)


#def adiabatic_expansion_due_to_mass_loss(a, M0, M1, m0, m1):
#    dmd = M1-M0
#    dma = m1-m0
#    print "\nStellar wind primary:", dmd, dma
#    # Note that mass loss rates are negative..
#    mdf = M1
#    maf = m1
#    mtf = mdf+maf
#    md = M0
#    ma = m0
#    mt = md+ma
#    dm  = mt-mtf
#    if dmd<zero:
#        alpha = dm/dmd
#        print "MT:", mt, mtf, md, mdf, ma, maf
##        fa = (((mdf/md) * (maf/ma))**(1/(1-alpha)) )**(-2) * (mt/mtf)
#        eta =  dma/dmd
#        a_new = a * pow(pow(mdf/md, eta)*maf/ma, -2)*mt/mtf
#        print "New a:", a, a_new
#        a = a_new
#        print "I'm not sure this is the correct use of the wind angular momentum loss equation"
#        #This equation is for one star loses mass in a wind and the other star accreting a fraction
#        #If a star loses wind mass, and the companion does not accrete: a_new = a * mt/mtf
#    return a
#
#def orbital_evolution_due_to_stellar_wind_mass_loss(bs):
#    a = adiabatic_expansion_due_to_mass_loss(bs.semimajor_axis, 
#                                             bs.child1.previous_mass, bs.child1.mass,
#                                             bs.child2.previous_mass, bs.child2.mass)
#    bs.semimajor_axis = a

def detached(tr):
    # orbital evolution is being taken into account in secular_code        
    if REPORT_FUNCTION_NAMES:
        print 'Detached'
        
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
                semi_detached(bs, bs.child1, bs.child2, outer_binary)
            elif not bs.child1.is_donor and bs.child2.is_donor:
                semi_detached(bs, bs.child2, bs.child1, tr)
            else:
                detached(tr)
                                        
        elif bs.child2.is_binary:
            if REPORT_BINARY_EVOLUTION:
                print bs.mass, bs.child1.mass, bs.child2.mass

            if bs.child1.is_donor:
                triple_mass_transfer()
            else:
                detached(tr)
                
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
                        
