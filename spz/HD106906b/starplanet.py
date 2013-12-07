"""
   example code for bridging a gravity solver with a hydrodynamics solver
"""
import numpy
from amuse.lab import *
from amuse.couple import bridge
from amuse import datamodel
from amuse.ext.evrard_test import uniform_unit_sphere

def get_kepler_elements(model_time, bh, star, converter):
    kep = Kepler(converter)
    kep.initialize_code()
    pos = bh.position - star.position
    vel = bh.velocity - star.velocity
    print "Kep:", bh.mass + star.mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]
    kep.initialize_from_dyn(bh.mass + star.mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    kep.stop()
    return a, e

def orbital_period(a, Mtot) :
    return 2*numpy.pi*(a**3/(constants.G*Mtot)).sqrt()

def construct_HD106906b(a, ecc):
    bodies = Particles(2)
    bodies[0].name = "HD106906"
    bodies[0].mass = 1.5| units.MSun
    bodies[0].radius = 1|units.RSun
    bodies[0].position = (0,0,0) | units.AU
    bodies[0].velocity = (0,0,0) | units.kms

    bodies[1].name = "HD106906b"
    bodies[1].mass = 11| units.MJupiter
    bodies[1].radius = 1|units.RJupiter
    bodies[1].position = (1,0,0) * a*(1+ecc)
    vc = (constants.G*bodies.mass.sum()/(a*(1+ecc))).sqrt()
    vc *= numpy.sqrt((1-ecc)/(1+ecc)) 

    bodies[1].velocity = (0,1,0) * vc
    bodies.move_to_center()
    return bodies

def gravity_hydro_bridge(a, ecc, t_end, n_steps):

    filename = "starlanet.hdf5"

    bodies = construct_HD106906b(a, ecc)
    print "Porb=", orbital_period(a, bodies.mass.sum())


    converter=nbody_system.nbody_to_si(bodies.mass.sum(), a)
    gravity = Huayno(converter)
    gravity.particles.add_particles(bodies)

    channel_from_gravity = gravity.particles.new_channel_to(bodies)
    channel_from_to_gravity = bodies.new_channel_to(gravity.particles)

    Ek0 = gravity.kinetic_energy
    Ep0 = gravity.potential_energy
    E0_tot = Ek0+Ep0

    dt = t_end/float(n_steps)
    time = zero
    while time < t_end:
        time += dt
        bodies.age = time

        write_set_to_file(bodies, filename, 'hdf5')

        gravity.evolve_model(time)
        channel_from_gravity.copy()
        a, ecc = get_kepler_elements(gravity.model_time, bodies[0], bodies[1], converter) 
        print "Time=", time, a, ecc
        Ek = gravity.kinetic_energy
        Ep = gravity.potential_energy
        E_tot = Ek+Ep
        print "Energies:", E_tot/E0_tot, Ek, Ep

    gravity.stop()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="int", default = 100,
                      help="number of diagnostics time steps [%default]")
    result.add_option("-a", unit=units.AU,
                      dest="a", type="float", default = 650|units.AU,
                      help="initial orbital separation [%default]")
    result.add_option("-e", dest="ecc", type="float", default = 0.9,
                      help="initial orbital eccentricity [%default]")
    result.add_option("-t", unit=units.yr, 
                      dest="t_end", type="float", default = 13000|units.yr,
                      help="end time of the simulation [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.AU, units.yr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")
    
    o, arguments  = new_option_parser().parse_args()
    gravity_hydro_bridge(**o.__dict__)
