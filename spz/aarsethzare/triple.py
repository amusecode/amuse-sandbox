from math import *
from numpy import *
import sys
import subprocess

import double_star
from EulerAngles import EulerAngles
from JacobiComCordinates import JacobiToCoMCoordinates
from massposvel_to_elements import *
from elements_to_massposvel import *
from amuse.legacy import *

def reorganize_units_quantitied(xa, unit) :
    xout = [] 
    print "Reorganize:", xa
    for i in range(len(xa)) :
        xout.append(xa[i].value_in(unit))
    xout = xout | unit
    print "xout:", xout
    return xout

def center_of_mass(mass, pos, vel) :
    print "Center of mass:", mass, pos, vel
    mtot = mass.sum()
    # reshape needed to prevent transpose the pos matrix.
    mr = mass.reshape((len(mass),1,))
    # sum() sums over the firs column.
    rcom = (mr*pos).sum(0)/mtot
    vcom = (mr*vel).sum(0)/mtot
    print "COM: ", rcom, vcom
    return rcom, vcom

def compress_array(x, y, z) :
    val = [[x[0].value_in(x.unit), y[0].value_in(x.unit), z[0].value_in(x.unit)], [x[1].value_in(x.unit), y[1].value_in(x.unit), z[1].value_in(x.unit)], [x[2].value_in(x.unit), y[2].value_in(x.unit), z[2].value_in(x.unit)]] |  x.unit
    return val

def extract_array(val) :
    x = [] | val.unit
    y = [] | val.unit
    z = [] | val.unit
    for i in range(len(val[0])) :
        x.append(val[i][0])
        y.append(val[i][1])
        z.append(val[i][2])
    return x, y, z

def aarseth_zare(t, mass, pos, vel, tend, unit_converter) :
    t = unit_converter.to_nbody(t)
    #mass = unit_converter.to_nbody(mass)
    pos = unit_converter.to_nbody(pos)
    vel = unit_converter.to_nbody(vel)
    tend = unit_converter.to_nbody(tend)
    print "aarseth_zare:", t, mass, pos, vel, tend
    ta = t.as_vector_with_length(3)
    tcrit = tend.as_vector_with_length(3)
    m = mass 
    x, y, z = extract_array(pos)
    vx, vy, vz = extract_array(vel)
    print "Before AZ", ta, m, x, y, z, vx, vy, vz, tcrit
    t, x, y, z, vx, vy, vz = instance.call_aarseth_zare(ta, m, x, y, z, vx, vy, vz, tcrit) 
    #t = tend
    print "After AZ", t, m, x, y, z, vx, vy, vz
    return t, x, y, z, vx, vy, vz

if __name__=="__main__":

    unit_mass = 1.0 | units.MSun
    unit_radius =  1. | units.RSun
    unit_converter = nbody_system.nbody_to_si(unit_mass, unit_radius)

    t = 0 | units.Myr
    if len(sys.argv) == 1 :
        m, pos, vel = generate_input()
    elif len(sys.argv) == 2 :
        filename = sys.argv[1]
        fptr = open(filename)
        m, pos, vel = read_massposvel(fptr)
        m = unit_converter.to_si(m)
        fptr.close()
    elif len(sys.argv) >=2 :
        print "command-line interface is not incorporated"
        exit()

    print "Print initial conditions (N-body units):", t, m, pos, vel

    bi, bo = massposvel_to_orbital_elements(m, pos, vel, unit_converter) 
    print "Binaries=", bi, bo
    Porb = bi.orbital_period()
    print "Orbital period", Porb, Porb.value_in(units.day)
    dmdt = 6E-10 * Porb.value_in(units.day) | units.MSun/units.Myr
    print "Mass transfer rate: ", dmdt

    Pinner = bi.orbital_period()
    Pouter = bo.orbital_period() 
#    dt = unit_converter.to_nbody(Pinner | units.day)
#    dt = unit_converter.to_nbody(Pinner | units.day)
#    dt = unit_converter.to_nbody(Pinner | units.day)
#    tend = unit_converter.to_si(dt).as_quantity_in(units.yr)
#    tend = 10.*Pinner | units.day
    tend = Pouter
    dt = Pinner
#    tend.as_quantity_in(units.yr)
    print "Orbital periods", dt, tend, Pinner, Pouter
#    dt = TP | nbody_system.time
    time = 0.0 | units.yr
    tinit = 0 | nbody_system.time
    while time < tend :
        print "mi=", bi.p.m, m[0], bi.p.m.value_in(units.MSun)
        m[0] = bi.p.m
        m[1] = bi.s.m
        m[2] = bo.s.m
        try :
            m = reorganize_units_quantitied(m, units.MSun)
        except :
            print "XM=", m
        print "M=", m
        m = unit_converter.to_nbody(m)
        print "M=", m
        print "Before AZ: Binaries=", bi, bo
        print tinit, m, pos, vel, dt
        tAZ, x, y, z, vx, vy, vz = aarseth_zare(tinit, m, pos, vel, dt, unit_converter)
        dt_yr = unit_converter.to_si(tAZ[0]).as_quantity_in(units.yr)
        time = time + dt_yr
        print "Time=", time, bi.t 

        pos = compress_array(x, y, z)
        vel = compress_array(vx, vy, vz)
        #bi, bo = construct_orbits(m, pos, vel)
        bi, bo = massposvel_to_orbital_elements(m, pos, vel, unit_converter)
        print "Binaries:", bi, bo

        bi.t = time.value_in(units.yr) 
        bo.t = time.value_in(units.yr) 
        print "Before MDot: Binaries=", bi, bo
        rcom, vcom = center_of_mass(m[:2], pos[:2], vel[:2])
        bi.rochelobe_overflow(bi.p, dt, dmdt)

        print "After MDot: Binaries=", bi, bo
        m, pos, vel = orbital_elements_to_massposvel(bi, bo)
        print "oe-M=", m

#        print_all_parameters(tAZ, x, y, z, vx, vy, vz)

    instance.stop()
