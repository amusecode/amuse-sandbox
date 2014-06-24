from amuse.lab import *
import numpy as np

REPORT_BINARY_EVOLUTION = False
REPORT_FUNCTION_NAMES = True

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



def semi_detached(bs, donor):
    if REPORT_FUNCTION_NAMES:
        print 'Semi detached'
    print "Roche-lobe overflow by star with mass:", donor.mass
    exit(0)
#    rstar = min(donor.radius, roche_radius(bs, donor))
#    Porb = orbital_period(bs)
#    mdc = donor.core_mass
#    Js = stellar_angular_momentum(donor, Porb)
#    Jb = binary_angular_momentum(bs)
#    print "Semi-detached:", Js, Jb, Jb/3.
#    if Js >= Jb/3. :
#        print "Darwin Riemann instability"
#        common_envelope(bs, donor)
#    else :
#        accretor = companion_star(bs, donor)
#        rochelobe_overflow(bs, donor, accretor)




def contact_binary():
    if REPORT_FUNCTION_NAMES:
        print "Contact binary"
    exit(0)


def adiabatic_expansion_due_to_mass_loss(a, M0, M1, m0, m1):
    dmd = M1-M0
    dma = m1-m0
    print "\nStellar wind primary:", dmd, dma
    # Note that mass loss rates are negative..
    mdf = M1
    maf = m1
    mtf = mdf+maf
    md = M0
    ma = m0
    mt = md+ma
    dm  = mt-mtf
    if dmd<zero:
        alpha = dm/dmd
        print "MT:", mt, mtf, md, mdf, ma, maf
#        fa = (((mdf/md) * (maf/ma))**(1/(1-alpha)) )**(-2) * (mt/mtf)
        eta =  dma/dmd
        a_new = a * pow(pow(mdf/md, eta)*maf/ma, -2)*mt/mtf
        print "New a:", a, a_new
        a = a_new
        print "I'm not sure this is the correct use of the wind angular momentum loss equation"
        #This equation is for one star loses mass in a wind and the other star accreting a fraction
        #If a star loses wind mass, and the companion does not accrete: a_new = a * mt/mtf
    return a

def orbital_evolution_due_to_stellar_wind_mass_loss(bs):
    a = adiabatic_expansion_due_to_mass_loss(bs.semimajor_axis, 
                                             bs.child1.old_mass, bs.child1.mass,
                                             bs.child2.old_mass, bs.child2.mass)
    bs.semimajor_axis = a

def detached(bs):
    if REPORT_FUNCTION_NAMES:
        print 'Detached'
#    pericenter = bs.semimajor_axis * (1-bs.eccentricity)
#    sum_of_radii = bs.child1.radius + bs.child2.radius
#    if bs.eccentricity>0 and pericenter < 3.*sum_of_radii :
#        tidal_circularization(bs)

    orbital_evolution_due_to_stellar_wind_mass_loss(bs)

def detached2(bs):
    #at the moment there are 2 functions called detached,
    # these functions do not take into account yet triple effects
    if REPORT_FUNCTION_NAMES:
        print 'Detached2'

#    pericenter = bs.semimajor_axis * (1-bs.eccentricity)
#    sum_of_radii = bs.child1.radius + bs.child2.radius
#    if bs.eccentricity>0 and pericenter < 3.*sum_of_radii :
#        tidal_circularization(bs)

    orbital_evolution_due_to_stellar_wind_mass_loss(bs)

def triple_mass_transfer():
    if REPORT_FUNCTION_NAMES:
        print 'Triple mass transfer'

def resolve_binary_interaction(bs):
   if bs.is_binary and bs.child1.is_star:
        Rl1 = roche_radius(bs, bs.child1)
        if REPORT_BINARY_EVOLUTION:
            print "Check for RLOF:", bs.child1.mass, bs.child1.old_mass
            print "Check for RLOF:", Rl1, bs.child1.radius
                
        if bs.is_binary and bs.child2.is_star:
            Rl2 = roche_radius(bs, bs.child2)
            if REPORT_BINARY_EVOLUTION:
                print "Check for RLOF:", bs.child2.mass, bs.child2.old_mass
                print "Check for RLOF:", Rl2, bs.child2.radius
            if Rl1 < bs.child1.radius: 
                if Rl2 < bs.child2.radius: 
                    contact_binary()
                else :
                    semi_detached(bs, bs.child1)
            else:
                if Rl2 < bs.child2.radius: 
                    semi_detached(bs, bs.child2)
                else:
                    detached(bs)
                                        
        elif bs.is_binary and bs.child2.is_binary:
            if REPORT_BINARY_EVOLUTION:
                print bs.mass, bs.child1.mass, bs.child2.mass
            if Rl1 < bs.child1.radius: 
                triple_mass_transfer()
            else:
                detached2(bs)
                
        else:
            print 'resolve binary interaction: type of system unknown'
            print bs.is_binary, 
            print bs.child1.is_binary, bs.child1.is_star, 
            print bs.child2.is_binary, bs.child2.is_star
            exit(-1) 
                               
   else:
        print 'resolve binary interaction: type of system unknown'
        print bs.is_binary, 
        print bs.child1.is_binary, bs.child1.is_star
#        print bs.child2.is_binary, bs.child2.is_star
        exit(-1)                    
                        
