from amuse.lab import *
from amuse.units.quantities import as_vector_quantity
from amuse.couple import encounters
from amuse.units import quantities
import numpy
import logging

def new_smalln(converter):
    result = SmallN(converter)
    result.parameters.timestep_parameter = 0.1
    result.parameters.cm_index = 2001
    return result

def new_kepler(converter):
    kepler = Kepler(converter)
    kepler.initialize_code()
    return kepler

def new_binary_orbit(mass1, mass2, semi_major_axis,
                     eccentricity = 0, keyoffset = 1):
    total_mass = mass1 + mass2
    mass_fraction_particle_1 = mass1 / (total_mass)
    
#    binary = Particles(keys=range(keyoffset, keyoffset+2))
    binary = Particles(2)
    binary[0].mass = mass1
    binary[1].mass = mass2
    
    mu = constants.G * total_mass
    
    velocity_perihelion = numpy.sqrt( mu / semi_major_axis  * ((1.0 + eccentricity)/(1.0 - eccentricity)))
    radius_perihelion = semi_major_axis * (1.0 - eccentricity)
    print velocity_perihelion
    
    binary[0].position = ((1.0 - mass_fraction_particle_1) * radius_perihelion * [1.0,0.0,0.0])
    binary[1].position = -(mass_fraction_particle_1 * radius_perihelion * [1.0,0.0,0.0])
    
    binary[0].velocity = ((1.0 - mass_fraction_particle_1) * velocity_perihelion * [0.0,1.0,0.0])
    binary[1].velocity = -(mass_fraction_particle_1 * velocity_perihelion * [0.0,1.0,0.0])

    return binary

def make_secondaries(center_of_masses, Nbin):

    resulting_binaries = Particles()
    singles_in_binaries = Particles()
    binaries = center_of_masses.random_sample(Nbin)
    for bi in binaries:
        mp = bi.mass
        ms = mp*numpy.random.random()

        a = 0.001 | units.parsec
        e = 0.6
        nb = new_binary_orbit(mp, ms, a, e) 
        nb.position += bi.position
        nb.velocity += bi.velocity
        nb = singles_in_binaries.add_particles(nb)
        nb.radius = 0.01 * a 
        bi.radius = 3*a 
        print "mas=", nb.mass, nb.key
        binary_particle = bi.copy()
        binary_particle.child1 = nb[0]
        binary_particle.child2 = nb[1]
        binary_particle.semi_major_axis = a
        binary_particle.eccentricity = e
        resulting_binaries.add_particle(binary_particle)

    single_stars = center_of_masses-binaries
    print single_stars
    print resulting_binaries
    print singles_in_binaries
    return single_stars, resulting_binaries, singles_in_binaries

def construct_orbital_elements(bi, converter):
    kep = new_kepler(converter)
    kep.set_longitudinal_unit_vector(1.0,0.0, 0.0)
    kep.set_transverse_unit_vector(0.0, 1.0, 0)
    comp1 = bi.child1
    comp2 = bi.child2
    mass = (comp1.mass + comp2.mass)
    pos = (comp2.position - comp1.position)
    vel = (comp2.velocity - comp1.velocity)
    kep.initialize_from_dyn(mass, pos[0], pos[1], pos[2],
                            vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    kep.stop()
    return a, e

def kira(tend, N, R, Nbin):
    logging.basicConfig(level=logging.ERROR)
    #encounters.LOG_ENERGY.setLevel(logging.DEBUG)

    mass = new_salpeter_mass_distribution(N, mass_min=10|units.MSun)
    converter = nbody_system.nbody_to_si(mass.sum(), R)
#    code = ph4(converter,redirection="none")
#    code = Huayno(converter)
    code = Hermite(converter)
    stars = new_plummer_model(N, convert_nbody=converter)
    stars.mass = mass
    stars.radius = 0.01/len(stars) | R.unit
    stars.velocity *= 0.
    stars[0].mass = 100 | units.MSun
    stars[0].position *= 0

    single_stars, binary_stars, singles_in_binaries = make_secondaries(stars, Nbin)
    print binary_stars

#    stellar = BSE()
    stellar = SeBa()
    stellar.particles.add_particles(single_stars)
    stellar.particles.add_particles(singles_in_binaries)
    stellar.binaries.add_particles(binary_stars)
    channel_to_stars = stellar.particles.new_channel_to(stars)

    encounter_code = encounters.HandleEncounter(
        kepler_code =  new_kepler(converter),
        resolve_collision_code = new_smalln(converter),
        interaction_over_code = None,
        G=constants.G
        )
    multiples_code = encounters.Multiples(
        gravity_code = code,
        handle_encounter_code = encounter_code,
        G=constants.G
        )
    multiples_code.particles.add_particles((stars-binary_stars).copy())
    multiples_code.singles_in_binaries.add_particles(singles_in_binaries)
    multiples_code.binaries.add_particles(binary_stars)
    multiples_code.commit_particles()
    channel_from_stars_to_particles = stellar.particles.new_channel_to(multiples_code.particles)

    stopping_condition = multiples_code.stopping_conditions.binaries_change_detection
    stopping_condition.enable()

    t = quantities.linspace(0*tend, tend, 11)
    for ti in t:
        print "t, Energy=", ti, multiples_code.particles.mass.sum(), multiples_code.get_total_energy()
        multiples_code.evolve_model(ti)
        print "at t=", multiples_code.model_time, "N-multiples:", len(multiples_code.multiples)

        if stopping_condition.is_set():
            new_binaries = stopping_condition.particles(0)
            lost_binaries = stopping_condition.particles(1)
            changed_binaries = stopping_condition.particles(2)

            # remove lost binaries before adding new ones (same star may appear in both).
            for bi in lost_binaries:
                print "remove old binary:", bi.key
                stellar.binaries.remove_particle(bi)

            for bi in new_binaries:
                print "add new binary:", bi
                a, e = construct_orbital_elements(bi, converter)
                bi.semi_major_axis = a
                bi.eccentricity = e
                stellar.binaries.add_particle(bi)
                print "new binary parameters", a, e
                print bi

            for bi in changed_binaries:
                bs = bi.as_particle_in_set(stellar.binaries)
                a, e = construct_orbital_elements(bi, converter)
                bs.semi_major_axis = a
                bs.eccentricity = e
                print "Modivide binarya arameters", a, e
                print bs

        stellar.evolve_model(ti)
        channel_from_stars_to_particles.copy_attributes(["mass", "radius"])

        kep = new_kepler(converter)
        kep.set_longitudinal_unit_vector(1.0,0.0, 0.0)
        kep.set_transverse_unit_vector(0.0, 1.0, 0)
# THIS NEEDS TO BE CHECKED
        print "Number of binaries=", len(stellar.binaries)
        for bi in stellar.binaries:
            bs = bi.as_particle_in_set(multiples_code.binaries)
            total_mass = bi.child1.mass+bi.child2.mass 
            kep.initialize_from_elements(total_mass,
                                         bi.semi_major_axis, bi.eccentricity)
            rel_position = as_vector_quantity(kep.get_separation_vector())
            rel_velocity = as_vector_quantity(kep.get_velocity_vector())
            mu = bi.child1.mass / total_mass 
            bs.child1.position = mu * rel_position 
            bs.child2.position = -(1-mu) * rel_position 
            bs.child1.velocity = mu * rel_velocity
            bs.child2.velocity = -(1-mu) * rel_velocity

            print "semi_major_axis=", ti, bi.semi_major_axis, total_mass, bi.child1.mass, bi.child2.mass, bi.eccentricity
        kep.stop()

        print "Lagrangian radii:", multiples_code.all_singles.LagrangianRadii(converter)
        print "MC.particles", multiples_code.particles
        print "Lagrangian radii:", multiples_code.particles.LagrangianRadii(converter)
        print "t, Energy=", ti, multiples_code.get_total_energy()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-t", unit=units.Myr,
                      dest="tend",type="float",default=1.|units.Myr)
    result.add_option("-R", unit=units.parsec,
                      dest="R",type="float",default=1|units.parsec)
    result.add_option("-N", 
                      dest="N",type="float",default=100)
    result.add_option("--Nbin", 
                      dest="Nbin",type="int",default=10)
    result.add_option("--seed", 
                      dest="seed",type="int",default=-1)
    return result

if __name__ == "__main__":
    set_printing_strategy("custom", 
                      preferred_units = [units.MSun, units.parsec, units.Myr], 
                      precision = 4, prefix = "", 
                      separator = " [", suffix = "]")

    options, arguments  = new_option_parser().parse_args()
    if options.seed>=0:
        numpy.random.seed(options.seed)
        # This is only for random.sample, which apparently does not use numpy
        import random
        random.seed(options.seed)
    kira(options.tend, options.N, options.R, options.Nbin)
