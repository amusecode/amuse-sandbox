"""
   example code for bridging a gravity solver with a hydrodynamics solver
"""
import numpy
from amuse.lab import *
from amuse.couple import bridge
from amuse import datamodel
from amuse.ext.evrard_test import uniform_unit_sphere

from starplanet import get_kepler_elements, orbital_period, construct_HD106906b

def calculate_disk_mass(disk_particles, starpos, Rmin, Rmax):
    
    inner_particles = disk_particles.select(lambda r: (starpos-r).length()<Rmin,["position"])
    outer_particles = disk_particles.select(lambda r: (starpos-r).length()>=Rmax,["position"])

    particles_in_annulus = disk_particles - inner_particles - outer_particles
    
    print "in/out disk:", len(disk_particles), len(inner_particles), len(outer_particles)
    print "Nannulus:", len(particles_in_annulus), particles_in_annulus.mass.sum()

    return particles_in_annulus.mass.sum()


def starplanetwithdisk(peri, apo, t_end, n_steps, Rinner, Router, Mdisk, Ndisk):
    a = 0.5*(peri+apo)
    ecc = (apo-peri)/(apo+peri)
    print "initial orbital parameters:", a, ecc
    bodies = construct_HD106906b(a, ecc)
    Pplanet = orbital_period(a, bodies.mass.sum())

    converter=nbody_system.nbody_to_si(bodies.mass.sum(), a)
    gravity = Huayno(converter) #, redirection="none")
    gravity.particles.add_particles(bodies)
    gravity.parameters.epsilon_squared = (1|units.AU)**2

    channel_from_gravity = gravity.particles.new_channel_to(bodies)
    channel_from_to_gravity = bodies.new_channel_to(gravity.particles)

    dt = t_end/float(n_steps)

    from amuse.ext.protodisk import ProtoPlanetaryDisk
    hydro_converter=nbody_system.nbody_to_si(bodies.mass.sum(), Rinner)
    disk_particles = ProtoPlanetaryDisk(Ndisk, convert_nbody=hydro_converter, 
                                densitypower=1.5, 
                                Rmin=1.0,
                                Rmax=Router/Rinner,
                                q_out=1.0,
                                discfraction=Mdisk/bodies.mass.sum()).result
    print "Mdisk=", Mdisk, disk_particles.mass.sum()
    bodies.move_to_center()
    com = disk_particles.center_of_mass()
    disk_particles.position += bodies[0].position
    disk_particles.velocity += bodies[0].velocity
#    ism.u = 0 | units.m**2 * units.s**-2 
#    ism.h_smooth= 0.01*a

    Pinnerdiskedge = orbital_period(Rinner, bodies[0].mass)
    print "Periods:", Pinnerdiskedge, Pplanet, orbital_period(Router, bodies[0].mass)

    hydro = Fi(converter) #, redirection="none")
    hydro.parameters.timestep = Pinnerdiskedge/1024.
    hydro.parameters.use_hydro_flag=True
    hydro.parameters.radiation_flag=False
    hydro.parameters.self_gravity_flag=True
    hydro.parameters.integrate_entropy_flag=False
    hydro.parameters.gamma=1.
    hydro.parameters.isothermal_flag=True
    hydro.parameters.epsilon_squared = (1|units.AU)**2
    hydro.gas_particles.add_particles(disk_particles)
    Eh0_tot = hydro.kinetic_energy + hydro.potential_energy + hydro.thermal_energy
    hydro.parameters.periodic_box_size = 10*a

    channel_from_hydro = hydro.gas_particles.new_channel_to(disk_particles)
    channel_to_hydro = disk_particles.new_channel_to(hydro.gas_particles)

    moving_bodies = ParticlesSuperset([bodies, disk_particles])

    filename = "starplanetdisk.hdf5"

    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(gravity, (hydro,) )
    gravhydro.add_system(hydro, (gravity,) )
    gravhydro.timestep = min(dt, min(0.1*Pplanet, 10*hydro.parameters.timestep))


    Ek0 = gravity.kinetic_energy
    Ep0 = gravity.potential_energy
    E0_tot = Ek0+Ep0

    time = 0 | units.Myr
    while time < t_end:
        time += dt
        bodies.age = time

        write_set_to_file(moving_bodies, filename, 'hdf5')

        gravhydro.evolve_model(time)
        channel_from_gravity.copy()
        channel_from_hydro.copy()

        a, ecc = get_kepler_elements(gravity.model_time, bodies[0], bodies[1], converter) 

        Mdisk = calculate_disk_mass(disk_particles, bodies[0].position, Rinner, Router)
        print "Diskmass=", time, Mdisk
        print "Planetary orbit=", time, a, ecc
        Ek = gravity.kinetic_energy
        Ep = gravity.potential_energy
        E_tot = Ek+Ep
        print "Energies:", E_tot/E0_tot, Ek, Ep


    gravhydro.stop()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="int", default = 10,
                      help="number of diagnostics time steps [%default]")
    result.add_option("-N", dest="Ndisk", type="int", default = 1024,
                      help="number of particles in the disk[%default]")
    result.add_option("-M", unit=units.MJupiter,
                      dest="Mdisk", type="float", default = 3.86984415126e-05|units.MJupiter,
                      help="Mass of the disk [%default]")
    result.add_option("-r", unit=units.AU,
                      dest="Rinner", type="float", default = 50|units.AU,
                      help="inner edge of the disk [%default]")
    result.add_option("-R", unit=units.AU,
                      dest="Router", type="float", default = 120|units.AU,
                      help="outer edge of the disk [%default]")
    result.add_option("--peri", unit=units.AU,
                      dest="peri", type="float", default = 50|units.AU,
                      help="pericenter of planet orbit [%default]")
    result.add_option("--apo", unit=units.AU,
                      dest="apo", type="float", default = 650|units.AU,
                      help="apocenter of planet orbit [%default]")
    result.add_option("-t", unit=units.yr, 
                      dest="t_end", type="float", default = 100|units.yr,
                      help="end time of the simulation [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.AU, units.yr], 
                          precision = 11, prefix = "", 
                          separator = " [", suffix = "]")
    o, arguments  = new_option_parser().parse_args()
    starplanetwithdisk(**o.__dict__)

