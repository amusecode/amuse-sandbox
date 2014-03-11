import os
import shutil
import cPickle

from amuse.units import nbody_system, units, constants
from amuse.io import read_set_from_file
from amuse.couple.bridge import CalculateFieldForCodesUsingReinitialize

from amuse.community.fi.interface import Fi
from amuse.community.octgrav.interface import Octgrav
from amuse.community.gadget2.interface import Gadget2
from amuse.community.ph4.interface import ph4
from amuse.community.sse.interface import SSE
from amuse.community.fastkick.interface import FastKick

from gravity_hydro_stellar import GravityHydroStellar



def simulate_young_star_cluster():
    end_time = 10.0 | units.Myr
    time_step = 0.01 | units.Myr
    
    new_working_directory()
    
    stars, gas, converter = load_young_star_cluster("YSC_stars10_gas1k_")
    #stars, gas, converter = load_young_star_cluster("YSC_stars1000_gas1000k_")
    system = new_gravity_hydro_stellar(stars, gas, converter, time_step)
    system.store_system_state()
    system.evolve_model(end_time)
    system.stop()

def new_working_directory():
    i = 0
    current_directory = os.getcwd()
    output_directory = os.path.join(current_directory, "output")
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    while os.path.exists(os.path.join(output_directory, "run_{0:=03}".format(i))):
        i += 1
    new_directory = os.path.join(output_directory, "run_{0:=03}".format(i))
    os.mkdir(new_directory)
    print "Created new directory for output:", new_directory
    os.mkdir(os.path.join(new_directory, "plots"))
    os.mkdir(os.path.join(new_directory, "snapshots"))
    shutil.copy(__file__, new_directory)
    os.chdir(new_directory)

def load_young_star_cluster(filename_base):
    filename_base = os.path.join("..", "..", "initial_conditions", filename_base)
    with open(filename_base + "info.pkl", "rb") as infile:
        converter = cPickle.load(infile)[0]
    stars = read_set_from_file(filename_base + "stars.amuse", "amuse")
    gas = read_set_from_file(filename_base + "gas.amuse", "amuse")
    return stars, gas, converter

def new_gravity_hydro_stellar(stars, gas, converter, time_step):
    stellar = new_stellar(stars)
    gravity = new_gravity("ph4", stars, converter)
    hydro = new_hydro(gas, converter, time_step)
    
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
        mode="gpu", 
        redirection="file", 
        redirect_file="grav_code.log")
    gravity.parameters.epsilon_squared = converter.to_si(0.001|nbody_system.length)**2
    print "Gravitational softening for stars:", gravity.parameters.epsilon_squared.sqrt().as_quantity_in(units.parsec)
    gravity.particles.add_particles(stars)
    gravity.commit_particles()
    return gravity

def new_gravity(code, *args, **kwargs):
    if code == "ph4":
        return new_ph4_gravity(*args, **kwargs)
    else:
        raise Exception("Undefined gravity code option: {0}".format(code))


def new_hydro(gas, converter, time_step):
    hydro = Gadget2(
        converter,
        node_label="cartesius", 
        number_of_workers=6, 
        redirection="file",
        redirect_file="gas_code.log")
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
    stellar = SSE(node_label="local")
    stellar.particles.add_particles(stars)
    return stellar

def new_gravity_field_from(codes, converter, epsilon_squared):
    gravity_field_code = Octgrav(
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

#count how often this function is called
new_gravity_field_from.counter = 1

if __name__ == "__main__":
    simulate_young_star_cluster()
