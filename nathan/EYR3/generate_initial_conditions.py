import os
import numpy
import cPickle

from amuse.units import nbody_system, units
from amuse.io import write_set_to_file

from amuse.ic.kroupa import new_kroupa_mass_distribution
from amuse.ic.fractalcluster import new_fractal_cluster_model
from amuse.ic.plummer import new_plummer_model
from amuse.ext.spherical_model import new_gas_plummer_distribution
from amuse.ext.relax_sph import relax
from amuse.community.fi.interface import Fi
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fastkick.interface import FastKick


def to_initial_conditions_directory():
    initial_conditions_directory = os.path.join(os.getcwd(), "initial_conditions")
    if not os.path.exists(initial_conditions_directory):
        os.mkdir(initial_conditions_directory)
        print "Created new initial conditions directory for output:", initial_conditions_directory
    os.chdir(initial_conditions_directory)

def generate_initial_conditions(
        number_of_stars = 10000,
        number_of_gas_particles = 2*10**6,
        star_formation_efficiency = 0.1,
        virial_radius = 0.33 | units.parsec,
        virial_ratio = 1.0,
        use_fractal = False):
    
    numpy.random.seed(12345678)
    seed_fractal = 312357271
    
    masses = new_kroupa_mass_distribution(number_of_stars)
    total_stellar_mass = masses.sum()
    total_mass = total_stellar_mass / star_formation_efficiency
    converter = nbody_system.nbody_to_si(total_mass, virial_radius)
    if use_fractal:
        stars = new_fractal_cluster_model(number_of_stars, convert_nbody=converter, do_scale=False, fractal_dimension=1.6, random_seed=seed_fractal)
    else:
        stars = new_plummer_model(number_of_stars, convert_nbody=converter, do_scale=False)
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
    
    filename = "YSC_{0}_stars{1}_gas{2}k_".format("fractal" if use_fractal else "plummer",
        number_of_stars, number_of_gas_particles/1000)
    print "Writing initial conditions to", filename, "+ stars/gas.amuse"
    write_set_to_file(stars, filename+"stars.amuse", "amuse", append_to_file=False)
    write_set_to_file(gas, filename+"gas.amuse", "amuse", append_to_file=False)
    with open(filename+"info.pkl", "wb") as outfile:
        cPickle.dump([converter], outfile)
    return stars, gas, filename

def new_hydro(gas, dynamical_timescale, converter):
    if False:
        hydro = Gadget2(converter, number_of_workers=8, redirection="file", redirect_file="gadget.log")
        hydro.parameters.time_max = 3 * dynamical_timescale
        hydro.parameters.max_size_timestep = dynamical_timescale / 100
        hydro.parameters.time_limit_cpu = 1.0 | units.Gyr
    else:
        hydro = Fi(converter, mode='openmp', redirection="file", redirect_file="fi.log")
        hydro.parameters.timestep = dynamical_timescale / 100
        hydro.parameters.eps_is_h_flag = True
    return hydro
    
def relax_initial_conditions(stars, gas, filename):
    dynamical_timescale = gas.dynamical_timescale()
    converter = nbody_system.nbody_to_si(dynamical_timescale, 1|units.parsec)
    hydro = new_hydro(gas, dynamical_timescale, converter)
    
    gravity_field_code = FastKick(converter, mode="gpu", number_of_workers=2)
    gravity_field_code.parameters.epsilon_squared = (0.01 | units.parsec)**2
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
#    relax_initial_conditions(stars, gas, filename)

