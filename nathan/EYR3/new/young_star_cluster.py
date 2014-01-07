import os
import shutil
import numpy
import time
import random

from amuse.units import nbody_system, units, constants
from amuse.couple.bridge import CalculateFieldForCodesUsingReinitialize

from amuse.community.fi.interface import Fi
from amuse.community.octgrav.interface import Octgrav
from amuse.community.gadget2.interface import Gadget2
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.ph4.interface import ph4
from amuse.community.sse.interface import SSE

from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.plummer import new_plummer_model
from amuse.ic.fractalcluster import new_fractal_cluster_model
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from gravity_hydro_stellar import GravityHydroStellar



def simulate_young_star_cluster():
    end_time = 10.0 | units.Myr
#~    end_time = 0.03 | units.Myr
    time_step = 0.01 | units.Myr
    
    new_working_directory()
    numpy.random.seed(12345678)
    random.seed(1234)
    
    stars, gas, converter = new_young_star_cluster()
    system = new_gravity_hydro_stellar(stars, gas, converter, time_step)
    system.evolve_model(end_time)
    
#~    current_time = time_step
#~    while system.model_time < end_time:
#~        system.evolve_model(current_time)
#~        current_time += time_step

def new_working_directory():
    i = 0
    current_directory = os.getcwd()
    while os.path.exists(os.path.join(current_directory, "run_{0:=03}".format(i))):
        i += 1
    new_directory = os.path.join(current_directory, "run_{0:=03}".format(i))
    os.mkdir(new_directory)
    print "Created new directory for output:", new_directory
    os.mkdir(os.path.join(new_directory, "plots"))
    os.mkdir(os.path.join(new_directory, "snapshots"))
    shutil.copy(__file__, new_directory)
    os.chdir(new_directory)

def new_young_star_cluster(
        number_of_stars = 2,
        number_of_gas_particles = 10000,
        star_formation_efficiency = 0.1,
        virial_radius = 1.0 | units.parsec,
        virial_ratio = 0.5):
    
    masses = new_salpeter_mass_distribution(number_of_stars, mass_max=100|units.MSun)
    masses = [20, 30.0] | units.MSun
    total_stellar_mass = masses.sum()
    total_mass = total_stellar_mass / star_formation_efficiency
    converter = nbody_system.nbody_to_si(total_mass, virial_radius)
#~    stars = new_plummer_model(number_of_stars, convert_nbody=converter)
    stars = new_fractal_cluster_model(number_of_stars, convert_nbody=converter, do_scale=False, fractal_dimension=1.6)
    stars.mass = masses
    stars.move_to_center()
    # scaling positions to match virial_radius:
    stars.position *= virial_radius / stars.virial_radius()
    # scaling velocities to match virial_ratio:
    stars.velocity *= numpy.sqrt(virial_ratio * converter.to_si(0.5|nbody_system.energy) * star_formation_efficiency / stars.kinetic_energy())
    
    gas = new_plummer_gas_model(
        number_of_gas_particles, 
        convert_nbody=converter, 
        base_grid=body_centered_grid_unit_cube)
    gas.h_smooth = 0.0 | units.parsec
    gas.mass = (total_mass - total_stellar_mass) / number_of_gas_particles
    # scaling positions to match virial_radius:
    gas.position *= virial_radius / gas.virial_radius()
    # scaling internal_energy to virial:
    gas.u *= converter.to_si(0.25|nbody_system.energy) * (1.0 - star_formation_efficiency) / gas.thermal_energy()
    return stars, gas, converter

def new_gravity_hydro_stellar(stars, gas, converter, time_step):
    gravity = new_gravity("ph4", stars, converter)
    hydro = new_hydro(gas, converter, time_step)
    stellar = new_stellar(stars)
    
    epsilon_squared_bridge = (constants.G * stars.mass.max() * time_step**2)**(2.0/3.0)
    print "Gravitational softening for Bridge:", epsilon_squared_bridge.sqrt().as_quantity_in(units.parsec)
    gravity_to_hydro = new_gravity_field_from([gravity], converter, epsilon_squared_bridge)
    hydro_to_gravity = new_gravity_field_from([hydro], converter, epsilon_squared_bridge)
    
    system = GravityHydroStellar(
        gravity, hydro, stellar, 
        gravity_to_hydro, hydro_to_gravity,
        time_step, 2*time_step, gas[0].mass)
    return system

def new_ph4_gravity(stars, converter):
    gravity = ph4(
        converter,
        node_label="GPU", 
#~        mode="gpu", 
        redirection="file", 
        redirect_file="grav_code.log")
    gravity.parameters.epsilon_squared = converter.to_si(0.001|nbody_system.length)**2
    print "Gravitational softening for stars:", gravity.parameters.epsilon_squared.sqrt().as_quantity_in(units.parsec)
#~    gravity.parameters.timestep_parameter = ?
    gravity.particles.add_particles(stars)
    gravity.commit_particles()
    return gravity

def new_gravity(code, *args, **kwargs):
    if code == "ph4":
        return new_ph4_gravity(*args, **kwargs)
    else:
        raise Exception("Undefined gravity code option: {0}".format(code))


def new_hydro(gas, converter, time_step):
#~    converter_for_hydro = nbody_system.nbody_to_si(1.0|units.Myr, converter.to_si(1.0|nbody_system.length))
    hydro = Gadget2(
        converter,
        node_label="CPU", 
        number_of_workers=10, 
        redirection="file",
        redirect_file="gas_code.log")
#    hydro = Gadget2(
#        converter,
#        node_label="cartesius", 
#        number_of_workers=48, 
#        number_of_nodes=2
#        redirection="file",
#        redirect_file="gas_code.log")
    hydro.parameters.gas_epsilon = converter.to_si(0.01|nbody_system.length) # or eps_is_h_flag
    hydro.parameters.n_smooth = 64
    hydro.parameters.n_smooth_tol = 0.005 |units.none
    hydro.parameters.time_max = 32.0 | units.Myr
    hydro.parameters.max_size_timestep = time_step
    hydro.parameters.time_limit_cpu = 1.0 | units.yr
    hydro.gas_particles.add_particles(gas)
    hydro.commit_particles()
    eps = hydro.gas_particles.radius.as_quantity_in(units.parsec)
    print "Gravitational softening for gas (mean, min, max):", eps.mean(), eps.min(), eps.max()
    return hydro

def new_stellar(stars):
    stellar = SSE()
    stellar.particles.add_particles(stars)
    return stellar

def new_gravity_field_from(codes, converter, epsilon_squared):
#~    gravity_field_code = Octgrav(
    gravity_field_code = BHTree(
        converter,
        node_label="GPU", 
        redirection="file", 
        redirect_file="gravity_field_code{0}.log".format(new_gravity_field_from.counter))
    gravity_field_code.parameters.epsilon_squared = epsilon_squared
    gravity_field_code.parameters.opening_angle = 0.5
    
    gravity_field = CalculateFieldForCodesUsingReinitialize(
        gravity_field_code, 
        codes, 
        required_attributes=['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz']
    )
    new_gravity_field_from.counter += 1
    return gravity_field
new_gravity_field_from.counter = 1


if __name__ == "__main__":
    new_young_star_cluster()

