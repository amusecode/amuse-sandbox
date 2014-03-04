import os
import numpy
import cPickle

from amuse.units import nbody_system, units
from amuse.io import write_set_to_file

from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.fractalcluster import new_fractal_cluster_model
from amuse.ext.spherical_model import new_gas_plummer_distribution
from amuse.ext.relax_sph import relax
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fastkick.interface import FastKick


def to_initial_conditions_directory():
    initial_conditions_directory = os.path.join(os.getcwd(), "initial_conditions")
    if not os.path.exists(initial_conditions_directory):
        os.mkdir(initial_conditions_directory)
        print "Created new initial conditions directory for output:", initial_conditions_directory
    os.chdir(initial_conditions_directory)

def generate_initial_conditions(
        number_of_stars = 10,
        number_of_gas_particles = 1000,
        star_formation_efficiency = 0.1,
        virial_radius = 1.0 | units.parsec,
        virial_ratio = 0.5):
    
    numpy.random.seed(12345678)
    seed_fractal = 312357271
    
    masses = new_salpeter_mass_distribution(number_of_stars, mass_max=100|units.MSun)
    total_stellar_mass = masses.sum()
    total_mass = total_stellar_mass / star_formation_efficiency
    converter = nbody_system.nbody_to_si(total_mass, virial_radius)
    stars = new_fractal_cluster_model(number_of_stars, convert_nbody=converter, do_scale=False, fractal_dimension=1.6, random_seed=seed_fractal)
    stars.mass = masses
    stars.move_to_center()
    print "scaling positions to match virial_radius"
    stars.position *= virial_radius / stars.virial_radius()
    print "scaling velocities to match virial_ratio"
    stars.velocity *= numpy.sqrt(virial_ratio * converter.to_si(0.5|nbody_system.energy) * star_formation_efficiency / stars.kinetic_energy())
    
    
    print "new_gas_plummer_distribution"
    gas = new_gas_plummer_distribution(
        number_of_gas_particles, 
        total_mass = (total_mass - total_stellar_mass), 
        virial_radius = virial_radius, 
        type = "fcc")
    gas.h_smooth = 0.0 | units.parsec
    
    filename = "YSC_stars{0}_gas{1}k_".format(number_of_stars, number_of_gas_particles/1000)
    write_set_to_file(stars, filename+"stars.amuse", "amuse")
    write_set_to_file(gas, filename+"gas.amuse", "amuse")
    with open(filename+"info.pkl", "wb") as outfile:
        cPickle.dump([converter], outfile)
    return stars, gas, filename

def relax_initial_conditions(stars, gas, filename):
    dynamical_timescale = gas.dynamical_timescale()
    converter = nbody_system.nbody_to_si(dynamical_timescale, 1|units.parsec)
    hydro = Gadget2(converter, number_of_workers=2)
    hydro.parameters.time_max = 3 * dynamical_timescale
    hydro.parameters.max_size_timestep = dynamical_timescale / 100
    hydro.parameters.time_limit_cpu = 1.0 | units.Gyr
    
    gravity_field_code = FastKick(converter)
    gravity_field_code.particles.add_particles(stars)
    relaxed_gas = relax(gas, hydro, gravity_field=gravity_field_code, 
        monitor_func="energy",
        bridge_options=dict(verbose=True, use_threading=False))
    gravity_field_code.stop()
    hydro.stop()
    write_set_to_file(relaxed_gas, filename+"gas_relaxed.amuse", "amuse")

if __name__ == "__main__":
    to_initial_conditions_directory()
    stars, gas, filename = generate_initial_conditions()
    relax_initial_conditions(stars, gas, filename)

