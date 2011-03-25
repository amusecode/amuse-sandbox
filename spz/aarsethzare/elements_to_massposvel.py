from math import *
from numpy import *
import sys
import subprocess

# prevent units from being printed:
#from amuse.support.data import console
#from amuse.support.data.console import set_printing_strategy
#set_printing_strategy(console.NoUnitsPrintingStrategy) 

import double_star
#import sys 
#sys.path.append("/home/spz/Instantaneous/AMUSE/amuse/sandbox/spz/binev/")
from EulerAngles import EulerAngles
from JacobiComCordinates import JacobiToCoMCoordinates

from amuse.community import *

def read_binary_line(line) :
    sline = line.split()
    ain = float(sline[7]) | units.RSun
    ein = float(sline[9]) 
    m1 = float(sline[11]) | units.MSun
    m2 = float(sline[13]) | units.MSun
    incl = float(sline[15])
    aperi = float(sline[16])
    llon = float(sline[17])
    iphase = float(sline[18])
#    b = Binary(ain, ein, m1, m2)
    b = double_star.double_star(ain, ein, m1, m2)
    b.incl = incl
    b.aperi = aperi
    b.llon = llon
    b.iphase = iphase
    return b

# read inner binary parameters: a, e, m1, m2, i, win, omin, phi1
def read_binary_parameters(fptr) :
    line = fptr.readline()
    sline = line.split("Binary")
    print sline
    bi = read_binary_line(sline[1]) 
    bo = read_binary_line(sline[2])
    return bi, bo

# define the inner orbit
# here the units are somewhat screwed-up, due to the goneometry
# therefore we split-off the units from the incoming double_star
def define_inner_orbit(b) : 
    a = b.a.value_in(units.RSun)
    e = b.e
    m1 = b.p.m.value_in(units.MSun)
    m2 = b.s.m.value_in(units.MSun)
    incl = radians(b.incl)
    aperi = radians(b.aperi)
    llon = radians(b.llon)
    phi = radians(b.phase)

    q=m2/m1
    m12=m1+m2
    Rpin=a*(1-e)
    Omegapin=sqrt((1+q)*(1+e)/Rpin**3)
    Tin=2*pi*sqrt(a**3/(1+q))

    rr=a*(1-e**2)/(1+e*cos(phi))
    rd=sqrt(m12/(a*(1-e**2)))*e*sin(phi)
    phid=sqrt(m12*a*(1-e**2))/rr**2
    r = array([rr*cos(phi), rr*sin(phi), 0])
    v = array([rd*cos(phi)-rr*phid*sin(phi), rd*sin(phi)+rr*phid*cos(phi), 0])

#   Rotate orbit and spins out of plane of outer orbit
    EA = array([incl, aperi, llon])
    r = EulerAngles(EA, r)
    v = EulerAngles(EA, v)
    return r, v

#call define_outer_orbit(aout, eout, m1, m2, m3, ophase)
def define_outer_orbit(b) :
    ophase = radians(b.phase)
    mtot = b.total_mass().value_in(units.MSun) # mass of inner binary and tertiariy (outer) star
    a = b.a.value_in(units.RSun)
    e = b.e

    rr=a*(1-e**2)/(1+e*cos(ophase))
    rd=sqrt(mtot/(a*(1-e**2)))*e*sin(ophase)
    phid=sqrt(mtot*a*(1-e**2))/rr**2
    r = array([rr*cos(ophase), rr*sin(ophase), 0])
    v = array([rd*cos(ophase)-rr*phid*sin(ophase), rd*sin(ophase)+rr*phid*cos(ophase), 0])
    return r, v

# Here in particular the unit conversion turned out to be rather tricky.
def orbital_elements_to_massposvel(bi, bo) :
    rin, vin = define_inner_orbit(bi)
    rout, vout = define_outer_orbit(bo)
    m = []
    m.append(bi.p.m.value_in(units.MSun))
    m.append(bi.s.m.value_in(units.MSun))
    m.append(bo.s.m.value_in(units.MSun))
    pos, vel = JacobiToCoMCoordinates(m, rin, vin, rout, vout) 
    pos = pos | units.RSun
    vel = vel | units.RSun / units.yr
    return m, pos, vel

# Dirty but effective.
# The assumption here is that the input is in Solar units and degrees
def read_input(sline) :
    ain = float(sline[1]) | units.RSun
    ein = float(sline[2])
    m1 = float(sline[3]) | units.MSun
    m2 = float(sline[4]) | units.MSun
    incl = float(sline[5])
    aperi = float(sline[6])
    llon = float(sline[7])
    iphase = float(sline[8])
    bi = double_star.double_star(ain, ein, m1, m2)
    bi.incl = incl
    bi.aperi = aperi
    bi.llon = llon
    bi.phase = iphase
    aout = float(sline[9]) | units.RSun
    eout = float(sline[10])
    peri = aout * (1-eout)
    m3 = float(sline[11]) | units.MSun
    ophase = float(sline[12])
    bo = double_star.double_star(aout, eout, m1+m2, m3)
    return bi, bo

if __name__=="__main__":

    # Define the set of units
    unit_mass      = 1.0 | units.MSun
    unit_radius    = 1.0 | units.RSun
    unit_converter = nbody_system.nbody_to_si(unit_mass, unit_radius)

    filename = 'elements_to_massposvel.in'
    print sys.argv
    if len(sys.argv) <= 2 :
        if len(sys.argv) == 1 :
            print "Default adopts input file called:", filename
        else :
            filename = sys.argv[1] 
        fptr = open(filename)
        bi, bo = read_binary_parameters(fptr)
        fptr.close()
    elif len(sys.argv) >2 :
        bi, bo = read_input(sys.argv)
        print "Generate relative cartesian coordinates for : "
        print bi, bo

    m, pos, vel = orbital_elements_to_massposvel(bi, bo)
#    for i in range(len(m)) :
#        print m[i], pos[i], vel[i]

    # for compatibility with it's sister function masposvel_to_elements 
    # this looks a bit awkward.
    for i in range(3) :
        st = str(i+1)+" "+str(m[i])+" "+str(pos[i][0].value_in(units.RSun))+" "+str(pos[i,1].value_in(units.RSun))+" "+str(pos[i,2].value_in(units.RSun))+" "+str(vel[i,0].value_in(units.RSun/units.yr))+" "+str(vel[i,1].value_in(units.RSun/units.yr))+" "+str(vel[i,2].value_in(units.RSun/units.yr))
        print st

