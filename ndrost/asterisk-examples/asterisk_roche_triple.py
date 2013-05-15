import os
import os.path

from matplotlib import pyplot
from amuse.units import units, nbody_system
from amuse.datamodel import Particles
from amuse.io import read_set_from_file
from amuse.ext.particles_with_color import new_particles_with_blackbody_color
from amuse.community.seba.interface import SeBa
from amuse.community.asterisk.interface import Asterisk

def start_visualisation(gas, stars):
    converter = nbody_system.nbody_to_si(1.0 | units.AU, gas.mass.sum() + stars.mass.sum())
    visualisation = Asterisk(converter, redirection="none")
    visualisation.initialize_code()
#    visualisation.parameters.use_octree_for_gas = True
#    visualisation.parameters.use_star_shader = False
    
    visualisation.gas_particles.add_particles(gas)
    visualisation.star_particles.add_particles(stars)
    return visualisation

def derive_temperatures_for(particles):
    instance = SeBa()
    in_code = instance.particles.add_particles(particles)
    instance.evolve_model(70|units.Myr)
    result = in_code.temperature
    instance.stop()
    return result

def add_visualization_attributes(particles, temperature_stars=None):
    if temperature_stars is None:
        particles.type = 2
        particles.radius = particles.h_smooth
    else:
        if temperature_stars == []:
            temperature_stars.extend(derive_temperatures_for(particles))
        particles.temperature = temperature_stars
        particles.type = 1
        particles.radius = 2 | units.RSun
    particles.alpha = 1.0
    new_particles_with_blackbody_color(particles)

def read_snapshot(gas_file, stars_file, directory, temperature_stars=[]):
    gas = read_set_from_file(os.path.join(directory, gas_file), "amuse")
    add_visualization_attributes(gas)
    stars = read_set_from_file(os.path.join(directory, stars_file), "amuse")
    add_visualization_attributes(stars, temperature_stars)
    return gas, stars

def visualize_simulation(directory):
    files = os.listdir(directory)
    files.sort()
    visualisation = start_visualisation(*read_snapshot(files[1], files[0], directory))
    step = 20
    for i in range(16):
        gas, stars = read_snapshot(files[2*i*step+1], files[2*i*step], directory)
        gas.u /= 10
        gas.synchronize_to(visualisation.gas_particles)
        gas.copy_values_of_attributes_to(["x", "y", "z", "red", "green", "blue"], visualisation.gas_particles)
        stars.synchronize_to(visualisation.star_particles)
        stars.copy_values_of_attributes_to(["x", "y", "z"], visualisation.star_particles)
        print "store_view"
        visualisation.store_view(i*0.4*10*step|units.day)
        print "done"
    visualisation.stop()

if __name__ in ('__main__', '__plot__'):
    #snapshot_directory = "/data1/vriesn/amuse/trunk/sandbox/nathan/roche_triple/validate/validate_4b/run_014_super_RLOF/snapshots"
    snapshot_directory = "/home/niels/Projects/AMUSE/roche_triple"
    visualize_simulation(snapshot_directory)
