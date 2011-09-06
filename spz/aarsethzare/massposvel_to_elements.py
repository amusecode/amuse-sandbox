from math import *
from numpy import *
import sys
import double_star
from amuse.community import *

#from amuse.support.datamodel import console
#from amuse.support.datamodel.console import set_printing_strategy
#set_printing_strategy(console.NoUnitsPrintingStrategy) 

from interface import AarsethZare 
instance = AarsethZare()

#      id  mass[Msun]        x y z [RSun]        vx vy vz [km/s]
spi = "0   1.0               0 0 0               0 0 0"
ssi = "1   3.00343905235e-06 1.49500000e+08 0 0  0 29.80 0"
sso = "2   3.24799388613e-07 2.27939100e+08 0 0  0 24.077 0"

def test_case_solar_system() :
    t = 0 | nbody_system.time
    mass = [1.0] | units.MSun
    mass.append(3.00343905235e-06 | units.MSun)
    mass.append(3.24799388613e-07| units.MSun)
    pos = [[0, 0, 0], [149.5e6, 0, 0],[227.94e6, 0, 0]] | units.km 
    vel = [[0, 0, 0], [0, 29.8000, 0],[0, 24.077, 0.0]] | units.km/units.s 
    return t, mass, pos, vel

def get_inner_binary_elements(m,pos,vel) :
    r = pos[1]-pos[0]
    v = vel[1]-vel[0]
    print "GIBE:", m, r, v
    e1, e2 = instance.construct_orbital_elements(m, r, v)
    return e1, e2

def get_outer_binary_elements(mass,pos,vel) :
    m12 = mass[0]+mass[1] 
    m123 = m12 + mass[2]
    q123 = m123/m12
    q123 = q123.value_in(q123.unit)
    m0 = 0.0  | nbody_system.mass
    print "Q=", q123, pos 
    print "Q*pos=", q123*pos[2] 
    r = []
    r.append(q123*pos[2].value_in(nbody_system.length))
    r = r | nbody_system.length
    r = r[0]
    v = []
    v.append(q123*vel[2].value_in(nbody_system.speed))
    v = v | nbody_system.speed
    v = v[0]
    print "v=", v
#    v.append(q123*vel[2])
    print "m=", m12, mass[2], m0
    m = []
    m.append(m12.value_in(m12.unit))
    m.append(mass[2].value_in(mass[2].unit))
    m.append(m0.value_in(m0.unit))
    m = m | nbody_system.mass
    print "GIBI:", m, r, v
    e1, e2 = instance.construct_orbital_elements(m,r,v)
#    elements = [m[0], m[1], e1[0], e1[1], e1[2], e2[0], e2[1], e2[2]] 
    return e1, e2

def massposvel_to_orbital_elements(m, pos, vel, unit_converter) :
    m = unit_converter.to_nbody(m)
    pos = unit_converter.to_nbody(pos)
    vel = unit_converter.to_nbody(vel)
    e1i, e2i = get_inner_binary_elements(m,pos,vel)
    print "Elements inner:", e1i, e2i
    bi = construct_binary(unit_converter.to_si(m), e1i, e2i)

    print "elements inner=", bi
    e1o, e2o = get_outer_binary_elements(m,pos,vel)
    bo = construct_binary(unit_converter.to_si(m), e1o, e2o)
    print "elements outer=", bo
    return bi, bo

def construct_binary(m, e1, e2) :
    print "construct_binary", m, e1, e2
    a = e1[0].value_in(e1[0].unit) | units.RSun
    e = e1[1].value_in(e1[1].unit)
    b = double_star.double_star(a, e, m[0], m[1])
    print "Binary=", b
    b.incl = degrees(e1[2].value_in(e1[2].unit))
    b.aperi = degrees(e2[0].value_in(e2[0].unit))
    b.llon = degrees(e2[1].value_in(e2[1].unit))
    b.phase = degrees(e2[2].value_in(e2[2].unit))
    return b

def read_binary_line(line) :

    #line = fptr.readline()
    sline = line.split()
    id = int(sline[0])
    m = float(sline[1]) 
    print "m=", m
    x = []
    x.append(float(sline[2]))
    x.append(float(sline[3]))
    x.append(float(sline[4]))
    print "x=", x
    v = [] 
    v.append(float(sline[5]))
    v.append(float(sline[6]))
    v.append(float(sline[7]))
    v  = v  
    print "v=", v
    return (id, m, x, v)

def generate_input() :
    inner_primary = read_binary_line(spi)
    inner_secondary = read_binary_line(ssi)
    outer_star = read_binary_line(sso)
    print "Inner primary:", inner_primary
    print "Inner secondary:", inner_secondary
    print "Outer star:", outer_star
    m = []
    m.append(inner_primary[1])
    m.append(inner_secondary[1])
    m.append(outer_star[1])
    m = m | units.MSun
    pos = [inner_primary[2], inner_secondary[2], outer_star[2]] | units.km
    vel = [inner_primary[3], inner_secondary[3], outer_star[3]] | units.km/units.s
    return m, pos, vel

def read_massposvel(fptr) :
    #line = fptr.readline()
    line = fptr.readline()
    inner_primary = read_binary_line(line)
    line = fptr.readline()
    inner_secondary = read_binary_line(line)
    line = fptr.readline()
    outer_star = read_binary_line(line)
    m = []
    m.append(inner_primary[1])
    m.append(inner_secondary[1])
    m.append(outer_star[1])
    m = m | nbody_system.mass
    pos = [inner_primary[2], inner_secondary[2], outer_star[2]] | nbody_system.length
    vel = [inner_primary[3], inner_secondary[3], outer_star[3]] | nbody_system.speed
    return m, pos, vel

if __name__=="__main__":

    # Define the set of units
    unit_mass = 1.0 | units.MSun
    unit_radius =  1. | units.RSun
    unit_converter = nbody_system.nbody_to_si(unit_mass, unit_radius)

    if len(sys.argv) == 1 :
        m, pos, vel = generate_input()
    elif len(sys.argv) == 2 :
        filename = sys.argv[1]
        fptr = open(filename)
        m, pos, vel = read_massposvel(fptr)
        fptr.close()
    elif len(sys.argv) >=2 :
        print "command-line interface is not incorporated"
        exit()

    print "MassPosVel:", m, pos, vel

    bi, bo = massposvel_to_orbital_elements(m, pos, vel, unit_converter) 
    print "Binaries", bi, bo
    print "Prob=", bi.orbital_period(), bo.orbital_period(), "Days"

    instance.stop()
