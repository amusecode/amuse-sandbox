"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
import numpy
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser


from amuse.community.kepler.interface import Kepler
def get_component_binary_elements(comp1, comp2, conv):
    kep = Kepler(redirection = "none")
    kep.initialize_code()

    mass = conv.to_nbody(comp1.mass + comp2.mass)
    pos = conv.to_nbody(comp2.position - comp1.position)
    vel = conv.to_nbody(comp2.velocity - comp1.velocity)
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    r = kep.get_separation()
    E,J = kep.get_integrals()	# per unit reduced mass, note
    kep.stop()

    return mass,a,e,r,E

def calculate_orbital_elements(gravity, converter):
    m,a,e,r,E = get_component_binary_elements(gravity.particles[0], gravity.particles[1], converter)
    m = converter.to_si(m).as_quantity_in(units.MSun)
    a = converter.to_si(a).as_quantity_in(units.AU)
    r = converter.to_si(r).as_quantity_in(units.AU)
    E = converter.to_si(r).as_quantity_in(units.AU)

    m0 = gravity.particles[0].mass
    m1 = gravity.particles[1].mass
    r_com = gravity.particles[0].position-gravity.particles[0].position
    v_com = gravity.particles[0].velocity-gravity.particles[0].velocity
    time = gravity.get_time().as_quantity_in(units.Myr)
    print "Orbital elements=", time, a, e, r, E, m, m0, m1, r_com, v_com
    return a, e

def orbital_period(a, Mtot) :
    return 2*numpy.pi*(a**3/(constants.G*Mtot)).sqrt()
    
def stellar_remnant_state(stellar):
    remnant = True
    for si in stellar.particles:
        if not (10 <= si.stellar_type.value_in(units.stellar_type) and \
                    si.stellar_type.value_in(units.stellar_type) < 16):
            remnant = False
    return remnant

def main(Mprim=3, Msec=2, a=100, e=0, t_end=10, nsteps=10, filename="wideWDb.hdf5", z=0.02, t_start=0, dtse_fraction=0.125):
    t_start = t_start | units.Myr
    t_end = t_end | units.Myr
    dt_diag = t_end/float(nsteps)
    dt_diag = -t_end
    Mprim = Mprim | units.MSun
    Msec = Msec | units.MSun
    a = a | units.AU

    stars=Particles(2)
    stars[0].mass= Mprim
    stars[1].mass= Msec
    Mstars = stars.mass.sum()
    print "Mstars=", stars.mass, Mstars

    stars[0].radius=0 | units.RSun
    stars[1].radius=0 | units.RSun

    stars[0].position = [0, 0, 0] | units.RSun
    stars[0].velocity = [0, 0, 0] | units.kms

    r = a
    v_orb = (constants.G*stars.mass.sum()*(2./r - 1./a)).sqrt().value_in(units.kms)

#    RLOF donor
    stars[1].position = [a.value_in(units.AU), 0, 0] | units.AU
    stars[1].velocity = [0, v_orb, 0] | units.kms

    converter=nbody_system.nbody_to_si(stars.mass.sum(),2*a)
#    converter=nbody_system.nbody_to_si(stars[0].mass,2*a)

    # star-up stellar evolution code
    stellar = SSE()
#    stellar = MESA()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(stars)
    stellar.commit_particles()

    stellar.evolve_model(t_start)

    channel_from_se_to_framework = stellar.particles.new_channel_to(stars)
    channel_from_se_to_framework.copy_attributes(["mass"])
    stars.move_to_center() 
#    print "stars: ", stars

    stars.scale_to_standard(convert_nbody=converter)

#    gravity = Hermite(converter)
#    gravity = ph4(converter)
    gravity = Huayno(converter)
    gravity.parameters.timestep_parameter = 0.01
    print "stars=", stars
    gravity.particles.add_particles(stars)

    channel_from_framework_to_gd = stars.new_channel_to(gravity.particles)
    channel_from_gd_to_framework = gravity.particles.new_channel_to(stars)
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    gravity.particles.move_to_center() 

    stars.a, stars.e = calculate_orbital_elements(gravity, converter)
    write_set_to_file(stars.savepoint(0|units.Myr), filename, 'hdf5')

    time = 0.0 | t_end.unit
#    P_orb = orbital_period(a, Mstars)
#    print "Orbital period=", P_orb
#    dt_se = dtse_fraction*P_orb
    t_diag = dt_diag
    while time<t_end:
    
        dt = 0.1*stellar.particles.time_step.amin()
        stellar.evolve_model(time+dt)
        time = stellar.model_time - t_start

        gravity.evolve_model(time)
        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

#        dt_gd = gravity.parameters.time_step
#        dt_gd = gravity.timestep
#        dt = min(dt_gd, dt_se)
#        print "timestep=", dt_gd, dt_se, dt

        channel_from_gd_to_framework.copy()
        channel_from_se_to_framework.copy_attributes(["mass", "radius"])
        channel_from_framework_to_gd.copy_attributes(["mass"])

        if time>t_diag:
            t_diag = time + dt_diag

            write_set_to_file(stars.savepoint(time), filename, 'hdf5')

            Ekin = gravity.kinetic_energy 
            Epot = gravity.potential_energy
            Etot = Ekin + Epot
            dE = Etot_prev-Etot
            dE_se = Etot_prev_se-Etot
            Mtot = stars.mass.sum()
            print "T=", time, 
            print "M=", Mtot, "(dM[SE]=", Mtot/Mstars, ")",
            print "E= ", Etot, "Q= ", Ekin/Epot,
            print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot, 
            print "(dE[SE]=", dE_se/Etot, ")"
#            print "stars: ", stars
            Etot_init -= dE
            Etot_prev = Etot
            stars.a, stars.e = calculate_orbital_elements(gravity, converter)
            peri = stars[-1].a*(1-stars[-1].e)
            
            if stellar_remnant_state(stellar):
                print "Both stars have turned into a remnant."
                exit()
            if peri <= stars.radius.sum():
                print "Stars touch at pericenter."
                exit()

    gravity.stop()
    stellar.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-n", dest="nsteps", type="float", default = 10,
                      help="diagnostics time steps [10]")
    result.add_option("-f", dest="filename", default = "widewdb.hdf5",
                      help="output filename [nbody.hdf5]")
    result.add_option("-M", dest="Mprim", type="float",default = 6,
                      help="Primary mass [3] MSun")
    result.add_option("--dtf", dest="dtse_fraction", type="float",default = 0.125,
                      help="se-gd time step fraction [0.125] Porb")
    result.add_option("-m", dest="Msec", type="float",default = 4,
                      help="secondary mass [2] MSun")
    result.add_option("-a", dest="a", type="float",default = 10000,
                      help="orbital separation [10000] in AU")
    result.add_option("-t", dest="t_end", type="float", default = 1.0,
                      help="end time of the simulation [0] Myr")
    result.add_option("--tstart", dest="t_start", type="float", default = 0.0,
                      help="start time for evolution [0] Myr")
    result.add_option("-e", dest="e", type="float", default = 0.0,
                      help="orbital eccentricity [0]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

