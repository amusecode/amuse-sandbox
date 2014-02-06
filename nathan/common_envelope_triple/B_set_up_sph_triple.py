import os
import os.path
import shutil
import numpy

from amuse.units import units, constants
from amuse.datamodel import Particles, Particle
from amuse.io import write_set_to_file
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH


def to_output_directory(directory_name="new"):
    model_directory = os.path.join(os.getcwd(), "initial_conditions", directory_name)
    if not os.path.exists(model_directory):
        if not os.path.exists(os.path.dirname(model_directory)):
            os.mkdir(os.path.dirname(model_directory))
        os.mkdir(model_directory)
        os.mkdir(os.path.join(model_directory, "plots"))
    shutil.copy(__file__, model_directory)
    os.chdir(model_directory)

def get_relative_velocity(total_mass, semimajor_axis, ecc):
    return (constants.G * total_mass * ((1.0 - ecc)/(1.0 + ecc)) / semimajor_axis).sqrt()

def triple_stars():
    triple = set_up_inner_binary()
    triple.add_particle(set_up_outer_star(triple.total_mass()))
    triple.move_to_center()
    return triple

def set_up_inner_binary():
    orbital_period = 835.0 | units.day
    masses = [10.0, 1.1] | units.MSun
    semimajor_axis = ( ((constants.G * masses.sum() * orbital_period**2) / 
        (4 * numpy.pi**2))**(1.0/3.0) ).as_quantity_in(units.RSun)
    
    print "   Initializing inner binary"
    print "   Semimajor axis inner binary:", semimajor_axis
    stars =  Particles(2)
    stars.mass = masses
    stars.position = [0.0, 0.0, 0.0] | units.AU
    stars.velocity = [0.0, 0.0, 0.0] | units.km / units.s
    stars[0].y = semimajor_axis
    stars[0].vx = -get_relative_velocity(stars.total_mass(), semimajor_axis, 0)
    stars.move_to_center()
    return stars

def set_up_outer_star(inner_binary_mass):
    orbital_period = 4020.0 | units.day
    mass = 1.3|units.MSun
    semimajor_axis = ( ((constants.G * (inner_binary_mass+mass) * orbital_period**2) / 
        (4 * numpy.pi**2))**(1.0/3.0) ).as_quantity_in(units.RSun)
    
    print "   Initializing outer star"
    print "   Semimajor axis outer binary:", semimajor_axis
    giant = Particle()
    giant.mass = mass
    giant.position = semimajor_axis * ([1, 0, 0] | units.none)
    giant.velocity = get_relative_velocity(giant.mass + inner_binary_mass, 
        semimajor_axis, 0) * ([0, 1, 0] | units.none)
    return giant

def giant_to_sph(giant, stellar_structure_file, number_of_sph_particles):
    model = convert_stellar_model_to_SPH(
        None, 
        number_of_sph_particles, 
        pickle_file = stellar_structure_file, 
        with_core_particle = True,
        target_core_mass  = 2.0 | units.MSun,
        do_store_composition = False,
        base_grid_options=dict(type="fcc")
    )
    sph_giant = model.gas_particles
    if model.core_particle is None:
        core = Particle(mass=0|units.MSun, radius=1|units.RSun)
        core_radius = core.radius
    else:
        core = model.core_particle
        core_radius = model.core_radius
    sph_giant.position += giant.position
    sph_giant.velocity += giant.velocity
    core.position = giant.position
    core.velocity = giant.velocity
    print "Core radius:", core_radius.as_string_in(units.RSun)
    return sph_giant, core, core_radius

def star_to_sph(star, number_of_sph_particles):
    stellar_evolution = stellar_evolution_code()
    se_star = stellar_evolution.particles.add_particle(star)
    stellar_evolution.evolve_model(22.7|units.Myr)
    model = convert_stellar_model_to_SPH(
        se_star, 
        number_of_sph_particles, 
        do_store_composition = False,
        base_grid_options=dict(type="fcc")
    )
    stellar_evolution.stop()
    sph_star = model.gas_particles
    sph_star.position += star.position
    sph_star.velocity += star.velocity
    return sph_star

    
def split_number_of_particles(number_of_sph_particles, triple):
    particle_mass = triple.total_mass() / number_of_sph_particles
    N2 = int(triple[1].mass / particle_mass + 0.5)
    N3 = int(triple[2].mass / particle_mass + 0.5)
    N1 = number_of_sph_particles - N2 - N3
    return N1, N2, N3


if __name__ == "__main__":
    number_of_sph_particles = 200000
    stellar_structure_file = os.path.join(os.getcwd(), "giant_models_MESA", "model_0008_836.7.pkl")
    
    to_output_directory("200k_837RSun")
    
    print "Initializing triple"
    triple = triple_stars()
    print "\nInitialization done:\n", triple
    
    N1, N2, N3 = split_number_of_particles(number_of_sph_particles, triple)
    
    print "Converting giant from file {0} to {1} SPH particles".format(stellar_structure_file, N1)
    sph_particles, core, core_radius = giant_to_sph(triple[0], stellar_structure_file, N1)
    
    print "Converting star to {0} SPH particles".format(N2)
    sph_particles.add_particles(star_to_sph(triple[1], N2))
    print "Converting star to {0} SPH particles".format(N3)
    sph_particles.add_particles(star_to_sph(triple[2], N3))
    print len(sph_particles)
    
    snapshotfile = "hydro_triple_gas.amuse"
    write_set_to_file(sph_particles, snapshotfile, format='amuse')
    snapshotfile = "hydro_triple_core.amuse"
    write_set_to_file(core.as_set(), snapshotfile, format='amuse')
