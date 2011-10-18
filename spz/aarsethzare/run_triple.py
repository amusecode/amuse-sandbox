from interface import AarsethZare 
#from numpy import *
from amuse.community import *

def get_outer_binary_elements(mass,pos,vel) :
    r = zeros(3)
    v = zeros(3)
    m12 = mass[0]+mass[1] 
    m123 = m12 + mass[2]
    q123 = m123/m12
    for i in range(3) :
        r[i]=q123*pos[i][2]
        v[i]=q123*vel[i][2]
    m = [m12, mass[2], 0.]
    e1, e2, error = instance.construct_orbital_elements(m,r,v)
    elements = [m[0],m[1], e1[0], e1[1], e1[2], e2[0], e2[1], e2[2]] 
    return elements

def get_inner_binary_elements(mass,pos,vel) :
    r = zeros(3)
    v = zeros(3)
    for i in range(3) :
        r[i]=pos[i][1]-pos[i][0]
        v[i]=vel[i][1]-vel[i][0]
    e1, e2, error = instance.construct_orbital_elements(mass,r,v)
    elements = [mass[0],mass[1], e1[0], e1[1], e1[2], e2[0], e2[1], e2[2]] 
    return elements

def aarseth_zare(particles, tend, unit_converter) :
    instance = AarsethZare(unit_converter, redirection="none")
    instance.particles.add_particles(particles)
    instance.evolve_model(tend)
    return instance.particles

def time_unit(elements) :
    mtot = elements[0]+elements[1]
    a = elements[2]
    import math
    #  pi=4.0d0*atan(1.0d0)
    orbital_time_unit = 2.*pi*sqrt(a**3./mtot)
    return orbital_time_unit

def DegenerateCoreMass(Porb, mmax=0.40, a=5.00, b=1.0e+5, c=0.110) :
    return min(mmax, (Porb/b)**(1./a) + c)

def check_orbits(m, pos, vel) :
    Einner = get_inner_binary_elements(m,pos,vel)
    print "elements inner=", Einner
    Eouter = get_outer_binary_elements(m,pos,vel)
    print "elements outer=", Eouter

    import sys
    sys.path.append("/home/spz/Latex/papers/2010/vdHvLNPZ/code/scr/")
    #from binary import *
    import binary

    bi = binary.Binary(Einner[2], Einner[3], Einner[0], Einner[1])
    bo = binary.Binary(Eouter[2], Eouter[3], Eouter[0], Eouter[1])
    tunit = time_unit(Eouter)
    year_per_step = bo.orbital_period()/tunit
    print "N-body time unit=", tunit
    return year_per_step

def SunEarthMars():
    particles = datamodel.Particles(3)
    particles.time = 0 | units.s
    particles.mass = [0, 5.9736e24, 6.46e23] | units.kg
    particles[0].mass = 1.0 | units.MSun
    particles.position = [[0,0,0],
                          [149.5e6, 0, 0],
                          [204.52e6, 0, 0]] | units.km 
    particles.velocity = [[0, 0, 0],
                          [0, 29800, 0],
                          [0, 24130, 0]] | units.m/units.s 
    return particles

if __name__=="__main__":

    unit_mass = 1 | units.MSun
    unit_radius = 1 | units.RSun
    unit_converter = nbody_system.nbody_to_si(unit_mass, unit_radius)

    particles = SunEarthMars()
    print particles

    tend = 100.0 | units.yr
    particles = aarseth_zare(particles, tend, unit_converter)
    print 't=', particles.time
    print "Particles:", particles

