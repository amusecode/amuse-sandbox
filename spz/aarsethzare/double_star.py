from random import *
from math import *
from amuse.community import *

class single_star :
    def __init__(self, m=0) :
        self.id = 0
        self.m = m
    def degenerate_coremass(self, Porb, mmax=0.40, a=5.00, b=1.0e+5, c=0.110) :
        return min(mmax, (Porb/b)**(1./a) + c) | units.MSun

class double_star :
    def __init__(self, a, e, mp, ms) :
        self.t = 0 | units.Myr
        self.a = a
        self.e = e
        self.p = single_star(mp)
        self.s = single_star(ms)
        self.incl = 0
        self.aperi = 0
        self.llon = 0
        self.phase = 0
    def __repr__(self) :
        tt = 'Binary (t= %s' % (self.t)
        tt += ' p= %s' % (self.p.id)
        tt += ' s= %s ' % (self.s.id)
        tt += ' (a= %s ' % (self.a)
        tt += 'e= %s ' % (self.e)
        tt += 'M= %s ' % (self.p.m)
        tt += 'm= %s ' % (self.s.m)
        tt += ' ( %s ' % (self.incl)
        tt += '%s ' % (self.aperi)
        tt += '%s ' % (self.llon)
        tt += '%s ))' % (self.phase)
        return tt
    def read_binary(self, line) :
        sline = line.split()
        print sline
        self.t = float(sline[2])
        self.p.id = float(sline[4])
        self.s.id = float(sline[6])
        self.a = float(sline[8])
        self.e = float(sline[10])
        self.p.m = float(sline[12])
        self.s.m = float(sline[14])
        self.incl = degrees(float(sline[16]))
        self.aperi = degrees(float(sline[17]))
        self.llon = degrees(float(sline[18]))
        self.phase = degrees(float(sline[19]))
    def orbital_period(self) :
        #print "OP:", self.a, self.p.m, self.s.m
        #a = self.a.value_in(units.RSun)
        #mtot = self.p.m.value_in(units.RSun) + self.s.m.value_in(units.RSun)
        #return 365.25 * sqrt((0.004652*a)**3/mtot) | units.day
        Pdays = 365.25 * sqrt((0.004652*self.a.value_in(units.RSun))**3/(self.p.m.value_in(units.MSun)+self.s.m.value_in(units.MSun))) #days
        return Pdays | units.day
    def roche_radius(self, p) :
        if p == self.p :
            return self.a * Rocheradius(self.p.m, self.s.m)
        else :
            return self.a * Rocheradius(self.s.m, self.p.m)
    def total_mass(self) :
        return self.p.m + self.s.m
    def companion_star(self, primary) :
        if primary == self.p :
            secondary = self.s
        else :
            secondary = self.p
        return secondary
    def rochelobe_overflow(self, donor, dt, dmdt) :
        accretor = self.companion_star(donor)
        print "Accretor:", accretor
        print "AMT:", dt, dmdt
        dmd =  dmdt * dt
        #dMEddington = 1.5e-8 * dt.value_in(dt.unit)
        MdotEddington = 1.5e-8 | units.MSun/units.yr
        dMEddington = MdotEddington * dt
        AccretionEfficiency = 1
        dma = min(AccretionEfficiency*dmd, dMEddington)
        print "Mdot: ", dmd, dma
        Pday = self.orbital_period().value_in(units.day)
        mcore = donor.degenerate_coremass(Pday)
        print "Donor:", mcore, "Pday=", Pday, donor.m, accretor.m, dma, dmd
        if donor.m-dmd < mcore :
            fm = max(0., (donor.m-mcore)/dmd)
            dmd = fm * dmd
            dma = fm * dma
            print 'Minimum mass reached for mass transfer', fm, dmd, dma
        self.transfer_mass(donor, dmd, dma)

    def transfer_mass(self, donor, dmd, dma) :
        accretor = self.companion_star(donor)
        print "Pre-RLOF:", donor.m, accretor.m, self.a
        md = donor.m
        ma = accretor.m
        mt = md + ma
        mdf = donor.m-dmd
        maf = accretor.m+dma
        mtf = mdf + maf
        dM  = mt-mtf
        alpha = (dM/dmd).value_in(units.none)
        fa = (((mdf/md) * (maf/ma))**(1/(1-alpha)) )**(-2) * (mt/mtf)
        print "fa=", fa
        self.a *= fa
        donor.m -= dmd
        accretor.m += dma
        print "Post-RLOF:", donor.m, accretor.m, self.a

def Rocheradius(M, m) :
    # Assure that the q is calculated in identical units.
    unit = M.unit
    # and that q itself has no unit
    q = M.value_in(unit)/m.value_in(unit)
    q13 =  q**(1./3.)
    q23 =  q13**2
    return  0.49*q23/(0.6*q23 + log(1 + q13))

def binding_energy(a, M, m):
  G = 6.67E-8
  Msun = 1.989E+33
  Rsun = 6.96e+10
  Eb = G * M *m * Msun**2/(2.*a*Rsun)
  return Eb

if __name__=="__main__":
    
    a = 1.0 | units.AU
    e = 0.000116700104236 | units.none
    M = 1.0 | units.MSun 
    m = 5.9736e24 | units.kg
    bs = double_star(a, e, M, m)

    print "double star =", bs
    print "Roche radius for the Sun:", bs.roche_radius(bs.p).as_quantity_in(units.AU)
    print "Roche radius for Earth:", bs.roche_radius(bs.s).as_quantity_in(units.AU)
    print "Orbital period: ", bs.orbital_period().as_quantity_in(units.yr)
