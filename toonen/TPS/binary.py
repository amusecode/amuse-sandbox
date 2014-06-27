from amuse.lab import *
import numpy as np
from math import sqrt

REPORT_BINARY_EVOLUTION = False
REPORT_FUNCTION_NAMES = False

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





def common_envelope_phase(bs, donor, accretor, tr):
    if REPORT_FUNCTION_NAMES:
        print 'Common envelope phase'

    #inner binary
    #remove mass
    #update star / se code
    #change inner orbit
    
    #outer binary
    #adiabatic_expansion_due_to_mass_loss ->instantaneous effect

    donor.is_donor = False
    # update stellar and orbital parameters in the secular code
    tr.channel_to_secular.copy()    
    

def stable_mass_transfer(bs, donor, accretor, tr):
    # orbital evolution is being taken into account in secular_code        
    if REPORT_FUNCTION_NAMES:
        print 'Stable mass transfer'

    #set mass transfer rate
    dt = tr.timestep
    dm = bs.mass_transfer_rate * dt
#    donor->self.se_code.particles[x].change_mass(dm, dt)
    
    
    Md = donor.mass
    Ma = accretor.mass
    print Md, Ma, 
    print donor.previous_mass, accretor.previous_mass
    # there is an implicit assumption in change_mass that the accreted mass is of solar composition (hydrogen)
#    accretor->self.se_code.particles[x].change_mass(-1*dm, dt)
    # for now, only conservative mass transfer   
#    accretor->self.se_code.particles[x].change_mass(-1*dm, -1*dt)
    triple.channel_from_se.copy()

    Md_new = donor.mass
    Ma_new = accretor.mass
    print Md_new, Ma_new, Md-Md_new, Ma-Ma_new
    print 'check if this is indeed conservative'
    accretion_efficiency = (Ma_new-Ma)/(Md-Md_new)
    print accretion_efficiency
#    bins[0].accretion_efficiency_mass_transfer = accretion_efficiency
    
    # update stellar parameters in the secular code
    print 'hmm where should info to secular code be updated?'
    tr.channel_to_secular.copy()   


def orbital_angular_momentum(bs):

    M = bs.child1.mass
    m = bs.child2.mass
    a = bs.semimajor_axis
    e = bs.eccentricity
    J = M*m * np.sqrt(constants.G*a*(1-e**2)/(M+m))

    if REPORT_BINARY_EVOLUTION:
        print 'Jorb:', M, m, a, e, J

    return J

def stellar_angular_momentum(ss, Porb):

    moment_of_inertia = ss.gyration_radius**2 * ss.mass * ss.radius**2
    Jstar = moment_of_inertia * ss.spin_angular_frequency
    return Jstar


def orbital_period(bs):

    Porb = 2*np.pi * np.sqrt(bs.semimajor_axis**3/constants.G / (bs.child1.mass + bs.child2.mass))
    return Porb


def mass_transfer_stability(bs, donor, accretor, tr):

    if REPORT_FUNCTION_NAMES:
        print 'Mass transfer stability'

    Porb = orbital_period(bs)
    Js = stellar_angular_momentum(donor, Porb)#accretor ook? check this
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
                
        if bs.is_binary and bs.child2.is_star:
            if REPORT_BINARY_EVOLUTION:
                Rl2 = roche_radius(bs, bs.child2)
                print "Check for RLOF:", bs.child2.mass, bs.child2.previous_mass
                print "Check for RLOF:", Rl2, bs.child2.radius

            if bs.child1.is_donor:            
                if bs.child2.is_donor:
                    contact_binary()
                else :
                    semi_detached(bs, bs.child1, bs.child2, tr)
            else:
                if bs.child2.is_donor: 
                    semi_detached(bs, bs.child2, bs.chil1, tr)
                else:
                    detached(tr)
                              
                                        
        elif bs.is_binary and bs.child2.is_binary:
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
                        
