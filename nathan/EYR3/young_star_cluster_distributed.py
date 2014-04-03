import os
import shutil
import cPickle
import logging
import cProfile

from amuse.units import nbody_system, units, constants
from amuse.io import read_set_from_file
from amuse.couple.bridge import CalculateFieldForCodesUsingReinitialize

from amuse.community.fi.interface import Fi
from amuse.community.octgrav.interface import Octgrav
from amuse.community.gadget2.interface import Gadget2
from amuse.community.ph4.interface import ph4
from amuse.community.sse.interface import SSE
from amuse.community.fastkick.interface import FastKick

from amuse.community.distributed.interface import DistributedAmuse, Resource, Pilot

from gravity_hydro_stellar import GravityHydroStellar



def simulate_young_star_cluster(distributed=None):
    end_time = 10.0 | units.Myr
    end_time = 0.06 | units.Myr
    time_step = 0.01 | units.Myr
    
    new_working_directory()
    
    stars, gas, converter = load_young_star_cluster("YSC_stars1000_gas100k_")
    #stars, gas, converter = load_young_star_cluster("YSC_stars1000_gas1000k_")
    system = new_gravity_hydro_stellar(stars, gas, converter, time_step, distributed)
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
#    gas = read_set_from_file(filename_base + "gas.amuse", "amuse")
    gas = read_set_from_file(filename_base + "gas_relaxed.amuse", "amuse")
    return stars, gas, converter

def new_gravity_hydro_stellar(stars, gas, converter, time_step, distributed):
    stellar = new_stellar(stars)
    gravity = new_gravity("ph4", stars, converter, distributed)
    hydro = new_hydro(gas, converter, time_step, distributed)
    
    epsilon_squared_bridge = (constants.G * stars.mass.max() * time_step**2)**(2.0/3.0)
    print "Gravitational softening for Bridge:", epsilon_squared_bridge.sqrt().as_quantity_in(units.parsec)
    gravity_to_hydro = new_gravity_field_from_stars([gravity], converter, epsilon_squared_bridge)
    hydro_to_gravity = new_gravity_field_from_gas([hydro], converter, epsilon_squared_bridge)
    
    system = GravityHydroStellar(
        gravity, hydro, stellar, 
        gravity_to_hydro, hydro_to_gravity,
        time_step, 2*time_step, gas[0].mass)
    return system

def new_ph4_gravity(stars, converter, distributed):
    if distributed is None:
        distributed_kwargs = dict()
    else:
        distributed_kwargs = dict(node_label="2GPUs_ph4", mode="gpu")
    
    gravity = ph4(
        converter,
        redirection="file", 
        redirect_file="grav_code.log",
        **distributed_kwargs)
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


def new_hydro(gas, converter, time_step, distributed):
    if distributed is None:
        distributed_kwargs = dict(number_of_workers=2)
    else:
        distributed_kwargs = dict(node_label="hydro", number_of_workers=24, number_of_nodes=1)
        if "cartesius" in distributed.resources.name:
            distributed.pilots.add_pilot(new_cartesius_pilot())
            print "Waiting for Cartesius reservation"
            distributed.wait_for_pilots()
    
    print "Start hydro"
    hydro = Gadget2(
        converter,
        redirection="file",
        redirect_file="gas_code.log",
        **distributed_kwargs)
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
    stellar = SSE(node_label="local", redirection="file", redirect_file="stellar.log")
    stellar.particles.add_particles(stars)
    return stellar

def new_gravity_field_from_stars(codes, converter, epsilon_squared):
    gravity_field_code = Octgrav(
        converter,
        node_label="GPU", 
        redirection="file", 
        redirect_file="gravity_field_code_gravity_to_hydro.log")
    gravity_field_code.parameters.epsilon_squared = epsilon_squared
    gravity_field_code.parameters.opening_angle = 0.5
    
    gravity_field = CalculateFieldForCodesUsingReinitialize(
        gravity_field_code, 
        codes, 
        required_attributes=['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'radius']
    ) # "radius" isn't really required, but get_state is faster than get_mass+get_position+get_velocity
    return gravity_field

def new_gravity_field_from_gas(codes, converter, epsilon_squared):
    gravity_field_code = Octgrav(
        converter,
        node_label="GPU", 
        redirection="file", 
        redirect_file="gravity_field_code_hydro_to_gravity.log")
    gravity_field_code.parameters.epsilon_squared = epsilon_squared
    gravity_field_code.parameters.opening_angle = 0.5
    
    gravity_field = CalculateFieldForCodesUsingReinitialize(
        gravity_field_code, 
        codes, 
        required_attributes=['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'u']
    ) # "u" isn't really required, but get_state_sph is faster than get_mass+get_position+get_velocity
    return gravity_field



def new_lgm_gateway():
    resource = Resource()
    resource.name = "LGM"
    resource.location = "fs.lgm.liacs.nl"
    resource.amuse_dir = "/home/vriesn/amuse-svn"
    return resource

def new_lgm_node(node_name):
    resource = Resource()
    resource.name = "LGM" + node_name
    resource.location = node_name
    resource.gateway = "LGM"
    resource.amuse_dir = "/home/vriesn/amuse-svn"
    return resource

def new_local_pilot():
    pilot = Pilot()
    pilot.resource_name = "local"
    pilot.node_count = 1
    pilot.time = 2 | units.hour
    pilot.slots_per_node = 10
    pilot.node_label = "local"
    return pilot

def new_gpu_node_pilot(resource, slots=2, label="GPU"):
    pilot = Pilot()
    pilot.resource_name = resource.name
    pilot.node_count = 1
    pilot.time = 2 | units.hour
    pilot.slots_per_node = slots
    pilot.node_label = label
    return pilot

def new_cpu_node_pilot(resource):
    pilot = Pilot()
    pilot.resource_name = resource.name
    pilot.node_count = 1
    pilot.time = 2 | units.hour
    pilot.slots_per_node = 10
    pilot.node_label = "CPU"
    return pilot

def new_cartesius_resource():
    resource = Resource()
    resource.name = "cartesius"
    resource.location = "int2-bb.cartesius.surfsara.nl"
    resource.amuse_dir = "/home/vriesn/amuse-svn"
    resource.scheduler_type = "slurm"
    return resource

def new_cartesius_pilot():
    pilot = Pilot()
    pilot.resource_name = "cartesius"
    pilot.node_count = 1
    pilot.time = 1 | units.hour
    pilot.queue_name = "short"
    pilot.slots_per_node = 24
    pilot.node_label = "hydro"
    return pilot

def new_hofvijver_resource():
    resource = Resource()
    resource.name = "hofvijver"
    resource.location = "hofvijver.strw.leidenuniv.nl"
    resource.amuse_dir = "/data1/vriesn/amuse-svn"
    return resource

def new_hofvijver_pilot():
    pilot = Pilot()
    pilot.resource_name = "hofvijver"
    pilot.node_count = 1
    pilot.time = 2 | units.hour
    pilot.slots_per_node = 64
    pilot.node_label = "hydro"
    return pilot


def start_distributed(lgm_node_names, forhydro="none"):
    lgm_nodes = [new_lgm_node(lgm_node_name) for lgm_node_name in lgm_node_names]
    instance = DistributedAmuse(redirection="file", redirect_file="distributed_amuse.log")
    instance.initialize_code()
    instance.resources.add_resource(new_lgm_gateway())
    for lgm_node in lgm_nodes:
        instance.resources.add_resource(lgm_node)
    
    for lgm_node in lgm_nodes:
        instance.pilots.add_pilot(new_cpu_node_pilot(lgm_node))
    instance.pilots.add_pilot(new_gpu_node_pilot(lgm_nodes[0], slots=1, label="2GPUs_ph4"))
    instance.pilots.add_pilot(new_gpu_node_pilot(lgm_nodes[1]))
    instance.pilots.add_pilot(new_local_pilot())
    
    if forhydro == "hofvijver":
        instance.resources.add_resource(new_hofvijver_resource())
        instance.pilots.add_pilot(new_hofvijver_pilot())
    elif forhydro == "cartesius":
        instance.resources.add_resource(new_cartesius_resource())
    else:
        raise Exception("Unknown option for forhydro resource: '{0}'".format(forhydro))
    
    print "Pilots:"
    print instance.pilots
    print "Waiting for pilots"
    instance.wait_for_pilots()
    return instance



if __name__ == "__main__":
#~    logging.basicConfig(filename='example.log', level=logging.WARN)
    logging.basicConfig(level=logging.WARN)
    logging.getLogger("code").setLevel(logging.DEBUG)

    run_local = False
    if run_local is True:
        cProfile.run("simulate_young_star_cluster()", "prof")
    else:
#~        distributed = start_distributed(lgm_node_names=["node05", "node07"], forhydro="cartesius")
        distributed = start_distributed(lgm_node_names=["node05", "node07"], forhydro="hofvijver")
        try:
            cProfile.run("simulate_young_star_cluster(distributed)", "prof")
        finally:
            distributed.stop()
            print "Done"

