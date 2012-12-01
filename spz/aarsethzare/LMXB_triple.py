from interface import AarsethZare 
instance = AarsethZare()
#from numpy import *
#Enable unit handling by import amuse.legacy
from amuse.lab import *
from math import *

import sys 
sys.path.append("/home/spz/Instantaneous/AMUSE/amuse/sandbox/spz/binev/")
import binary

def construct_binary(m, e1, e2) :
    m1 = m[0].value_in(m[2].unit) #| units.MSun
    m2 = m[1].value_in(m[2].unit) #| units.MSun
    a = e1[0].value_in(e1[0].unit) #| units.RSun
    e = e1[1].value_in(e1[1].unit)
    b = binary.Binary(a, e, m1, m2)
    b.incl = degrees(e1[2].value_in(e1[2].unit))
    b.aperi = degrees(e2[0].value_in(e2[0].unit))
    b.llon = degrees(e2[1].value_in(e2[1].unit))
    b.phase = degrees(e2[2].value_in(e2[2].unit))
    return b

def get_outer_binary_elements(mass,pos,vel) :
    m = [] | nbody_system.mass
    r = [] | nbody_system.length
    v = [] | nbody_system.speed
    m12 = mass[0]+mass[1] 
    m123 = m12 + mass[2]
    q123 = m123/m12
    m0 = 0.0  | nbody_system.mass
    #print "Q=", q123
    r.append(q123*pos[2])
    v.append(q123*vel[2])
    m.append(m12)
    m.append(mass[2])
    m.append(m0)
    print "GIBI:", m, r, v
    e1, e2 = instance.construct_orbital_elements(m,r,v)
#    elements = [m[0], m[1], e1[0], e1[1], e1[2], e2[0], e2[1], e2[2]] 
    bo = construct_binary(m, e1, e2)

    return bo

def get_inner_binary_elements(m,pos,vel) :
    r = [] | nbody_system.length
    v = [] | nbody_system.speed
    r.append(pos[1]-pos[0])
    v.append(vel[1]-vel[0])
    print "GIBE:", m, r, v
    e1, e2 = instance.construct_orbital_elements(m,r,v)
    bi = construct_binary(m, e1, e2)
#    elements = [m[0], m[1], e1[0].value_in(e1[0].unit) | units.RSun, e1[1], e1[2], e2[0], e2[1], e2[2]] 
    return bi

def test_case_solar_system() :
    t = 0 | nbody_system.time
    mass = [1.0] | units.MSun
    mass.append(5.9736e24 | units.kg)
    mass.append(6.46e23 | units.kg)
    pos = [[0, 0, 0], [149.5e6, 0, 0],[204.52e6, 0, 0]] | units.km 
    vel = [[0, 0, 0], [0, 29800.0, 0],[0, 24130., 0.0]] | units.m/units.s 
    return t, mass, pos, vel

def construct_initial_conditions() :
    #Inner: Binary (t= 0 p= 0 s= 0 (a= 52.8132361507 e= 0 M= 1.1 m=
    #1.3 ( 42.4911139342 221.337219619 25.0450782222 274.667014114 ))
    #Outer: Binary (t= 0 p= 0 s= 0 (a= 1400.97007777 e=
    #0.678139318564 M= 2.4 m= 1.0 ( 0 0 0 79.0583859092 ))
    t = 0 | nbody_system.time
    mass = [1.1, 1.3, 1.0] | nbody_system.mass #secretly MSun
    pos = [[-12.578787679, -58.4511049545, 89.8231028877], #xxx
           [-198.157193278, -189.734002337, 464.627115644],#yyy
           [-13.4221707428, 11.3572213977, 0.0]] | nbody_system.length#zzz
    vel = [[0.0660907760075, -0.0201875848322, -0.0464559933263],
           [0.0723351192845, -0.0927973717159, 0.0410679520178],
           [0.0561101265997, -0.0474777994305, 0.0]] | nbody_system.speed
    return t, mass, pos, vel

def extract_array(val) :
    x = [] | val.unit
    y = [] | val.unit
    z = [] | val.unit
    for i in range(len(val[0])) :
        x.append(val[i][0])
        y.append(val[i][1])
        z.append(val[i][2])
    return x, y, z

def split_off_units(x) :
    y = [] |  x[0].unit
    for i in x :
        y.append(i.value_in(x[0].unit))
    return y

def combine_arrays_to_matrix(x, y, z) :
    val = [[split_off_units(x), split_off_units(y), split_off_units(z)]]
    return val

def compress_array(x, y, z) :
    val = [[x[0].value_in(x.unit), y[0].value_in(x.unit), z[0].value_in(x.unit)], [x[1].value_in(x.unit), y[1].value_in(x.unit), z[1].value_in(x.unit)], [x[2].value_in(x.unit), y[2].value_in(x.unit), z[2].value_in(x.unit)]] |  x.unit
    return val

def aarseth_zare(t, mass, pos, vel, tend) :
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

def DegenerateCoreMass(Porb, mmax=0.40, a=5.00, b=1.0e+5, c=0.110) :
    return min(mmax, (Porb/b)**(1./a) + c)

def construct_orbits(m, pos, vel) :
    bi = get_inner_binary_elements(m,pos,vel)
    print "elements inner=", bi
    bo = get_outer_binary_elements(m,pos,vel)
    print "elements outer=", bo
    return bi, bo

def read_triple(infile) :
    id = []
    m = []
    pos = []
    vel = []
    fptr = open(infile)
    line = fptr.readline()
    t = float(line.split()[0]) | nbody_system.time
    line = fptr.readline()
    while line :
        sline = line.split()
        id.append(int(sline[0]))
        m.append(float(sline[1]))
        pos.append([float(sline[2]), float(sline[3]), float(sline[4])])
        vel.append([float(sline[5]), float(sline[6]), float(sline[7])])
        line = fptr.readline()
    fptr.close()
    m = m | nbody_system.mass
    print "Pre pos: ",  pos
    pos = pos | nbody_system.length
    print "Pos pos:", pos
    vel = vel | nbody_system.speed
    return t, m, pos, vel

def print_all_parameters(tAZ, x, y, z, vx, vy, vz) :
    print 't=', tAZ[0], 'r=', x, y, z, "v=", vx, vy, vz
    print 't(Units)=', unit_converter.to_si(tAZ[0]).as_quantity_in(units.yr)
    print "total mass = ", unit_converter.to_si(m).sum().as_quantity_in(units.MSun)
    print "r=", unit_converter.to_si((x**2+y**2).sqrt()).as_quantity_in(units.AU)
    print "v=", unit_converter.to_si((vx**2+vy**2).sqrt()).as_quantity_in(units.km/units.s)
    print "x=", unit_converter.to_si(x).as_quantity_in(units.AU)
    print "y=", unit_converter.to_si(y).as_quantity_in(units.AU)
    print "z=", unit_converter.to_si(z).as_quantity_in(units.AU)
    print "vx=", unit_converter.to_si(vx).as_quantity_in(units.km/units.s)
    print "vy=", unit_converter.to_si(vy).as_quantity_in(units.km/units.s)
    print "vz=", unit_converter.to_si(vz).as_quantity_in(units.km/units.s)

def ApplyMassTransfer(bi, dt, dmdt) :
    print "AMT:", dt.value_in(dt.unit), dmdt
    dmd =  dmdt * dt.value_in(dt.unit)
    dMEddington = 1.5e-8 * dt.value_in(dt.unit)
    AccretionEfficiency = 1
    dma = min(AccretionEfficiency*dmd, dMEddington)
    Porb = bi.orbital_period()
    mcore = DegenerateCoreMass(Porb)
    print mcore, Porb, bi.p.m, dmd
    if bi.p.m-dmd < mcore :
        fm = max(0., (bi.p.m-mcore)/dmd)
        dmd = fm * dmd
        dma = fm * dma
        print 'Minimum mass reached for mass transfer', fm, dmd, dma, mass
    bi.evolve_nonconservative(dmd, dma)
    return bi

def ncenter_of_mass(mass, pos, vel) :
    rcom = [] 
    vcom = [] 
    for i in range(3) :
        rcom.append(mass[0]*pos[0][i] + mass[1]*pos[1][i])
        vcom.append(mass[0]*vel[0][i] + mass[1]*vel[1][i])
    mtot = mass[0]+mass[1]
    for i in range(len(rcom)) :
        rcom[i] = rcom[i]/mtot
        vcom[i] = vcom[i]/mtot
#    rcom = rcom | nbody_system.length
#    vcom = vcom | nbody_system.speed
    return rcom, vcom

def mcenter_of_mass(mass, pos, vel) :
    print "Center of mass:", mass, pos, vel
    rcom = [] 
    vcom = [] 
    for i in range(3) :
        rcom.append(mass[0]*pos[0][i] + mass[1]*pos[1][i])
        vcom.append(mass[0]*vel[0][i] + mass[1]*vel[1][i])
    mtot = mass.sum()
    for i in range(3) :
        rcom[i] = rcom[i]/mtot
        vcom[i] = vcom[i]/mtot
    rcom = rcom | nbody_system.length
    vcom = vcom | nbody_system.speed
    print "Com:", rcom, vcom
    return rcom, vcom

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

def reconstruct_relative_position_and_velocity(bi, pos, vel, rcom, vcom) :

    m0 = 0.0 | nbody_system.mass
    m = [] | nbody_system.mass
    m.append(bi.p.m | nbody_system.mass)
    m.append(bi.s.m | nbody_system.mass)
    m.append(m0)
    print "XXXX", m, pos, vel, rcom, vcom
    e1 = [bi.a, bi.e, bi.incl]
    e2 = [bi.aperi, bi.llon, bi.phase]
    print "XXXX", e1, e2
    r, v = instance.construct_orbital_coordinates(m, e1, e2)
    print "Post Transformation: r, v", r, v
    mtot = m[0]+m[1]
    for i in range(3) :
        pos[0][i] = rcom[i] - r[i]* m[1]/mtot
        vel[0][i] = vcom[i] - v[i]* m[1]/mtot
        pos[1][i] = rcom[i] + r[i]* m[0]/mtot
        vel[1][i] = vcom[i] + v[i]* m[0]/mtot
    return pos, vel

def temporary_storage_function_not_used() :
        sys.path.append("/home/spz/Instantaneous/AMUSE/amuse/sandbox/spz/binev/")
        #from EulerAngles import EulerAngles
        from JacobiComCordinates import JacobiToCoMCoordinates
        from elementstoposvel import *

        rin, vin = define_inner_orbit((bi.a, bi.e, bi.p.m, bi.s.m, bi.incl, bi.aperi, bi.llon, bi.phase))
        mm = [bi.p.m, bi.s.m, bo.s.m]
        rout, vout = define_outer_orbit(bo.a, bo.e, mm, bo.phase) 
        pos, vel = JacobiToCoMCoordinates(mm, rin, vin, rout, vout) 
        print "Pre reconstruction: ", mm, pos, vel
        m = mm | nbody_system.mass
        pos = convert_to_proper_units(pos, nbody_system.length)
        vel = convert_to_proper_units(vel, nbody_system.length/nbody_system.time)
        print "Post reconstruction: ", m, pos, vel
        bi, bo = construct_orbits(m, pos, vel)
        print "the resulting binaries: ", bi, bo


def convert_to_proper_units(pos, unit) :
    p = [] 
    for i in range(len(pos[0])) :
        print "i=", i
        p.append([pos[i][0],pos[i][1],pos[i][2]])
    return p | unit

if __name__=="__main__":

    unit_mass = 1.0 | units.MSun
    unit_radius =  1. | units.RSun
    unit_converter = nbody_system.nbody_to_si(unit_mass, unit_radius)

    if len(sys.argv) > 1 :
        filename = sys.argv[1]
        t, m, pos, vel = read_triple(filename)
    else :
        t, m, pos, vel = test_case_solar_system()
        m = unit_converter.to_nbody(m)
        pos = unit_converter.to_nbody(pos)
        vel = unit_converter.to_nbody(vel)

    print "Print initial conditions (N-body units):", t, m, pos, vel

    bi, bo = construct_orbits(m, pos, vel)
    print "Binaries=", bi, bo
    Porb = bi.orbital_period()
    dmdt = 6E-10 * Porb
    print "Mass transfer rate: ", dmdt

    VP2 = 2.0*(bo.p.m + bo.s.m)/bo.a
    TP = bo.a/sqrt(VP2)

    rmax = max(bi.a,bo.a)
    TSTAR = rmax*sqrt(rmax/(bi.p.m + bi.s.m + bi.s.m))
    print "Aarseth TP=", TP, TSTAR
    TSTAR = min(TP,TSTAR)

    Pinner = bi.orbital_period()
    Pouter = bo.orbital_period() 
#    dt = unit_converter.to_nbody(Pinner | units.day)
#    dt = unit_converter.to_nbody(Pinner | units.day)
    dt = unit_converter.to_nbody(Pinner | units.day)
#    tend = unit_converter.to_si(dt).as_quantity_in(units.yr)
    tend = 10.*Pinner | units.day
    tend.as_quantity_in(units.yr)
    print "Orbital periods", dt, tend, Pinner, Pouter
#    dt = TP | nbody_system.time
    time = 0.0 | units.yr
    while time < tend :
        print "Before AZ: Binaries=", bi, bo
        tAZ, x, y, z, vx, vy, vz = aarseth_zare(t, m, pos, vel, dt)
        dt_yr = unit_converter.to_si(tAZ[0]).as_quantity_in(units.yr)
        time = time + dt_yr
        print "Time=", time, bi.t 

        pos = compress_array(x, y, z)
        vel = compress_array(vx, vy, vz)
        bi, bo = construct_orbits(m, pos, vel)
        bi.t = time.value_in(units.yr) 
        bo.t = time.value_in(units.yr) 
        print "Before MDot: Binaries=", bi, bo
        rcom, vcom = center_of_mass(m[:2], pos[:2], vel[:2])
        #bi = ApplyMassTransfer(bi, dt_yr, dmdt)
        print "After MDot: Binaries=", bi, bo

        pos, vel = reconstruct_relative_position_and_velocity(bi, pos, vel, rcom, vcom)
#        print_all_parameters(tAZ, x, y, z, vx, vy, vz)
        m[0] = bi.p.m | nbody_system.mass
        m[1] = bi.s.m | nbody_system.mass
        m[2] = bo.s.m | nbody_system.mass

    instance.stop()
