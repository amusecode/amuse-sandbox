"""
 run script for clustergas runs using DistributedAmuse

"""

import logging
import cProfile

from amuse.community.distributed.interface import DistributedAmuse, Resource, Pilot
from amuse.lab import *

from young_star_cluster import simulate_young_star_cluster

def new_lgm_gateway():
    resource = Resource()
    resource.name = "LGM"
    resource.location = "fs.lgm.liacs.nl"
    resource.amuse_dir = "/home/niels/amuse"
    return resource

def new_lgm_node(node_name):
    resource = Resource()
    resource.name = "LGM" + node_name
    resource.location = node_name
    resource.gateway = "LGM"
    resource.amuse_dir = "/home/niels/amuse"
    return resource

def new_local_pilot():
    pilot = Pilot()
    pilot.resource_name = "local"
    pilot.node_count = 1
    pilot.time = 2 | units.hour
    pilot.slots_per_node = 10
    pilot.node_label = "local"
    return pilot

def new_gpu_node_pilot(resource):
    pilot = Pilot()
    pilot.resource_name = resource.name
    pilot.node_count = 1
    pilot.time = 2 | units.hour
    pilot.slots_per_node = 2
    pilot.node_label = "GPU"
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
    resource.location = "ndrosta@int2-bb.cartesius.surfsara.nl"
    resource.amuse_dir = "/home/ndrosta/amuse-svn"
    resource.scheduler_type = "slurm"
    return resource

def new_cartesius_pilot():
    pilot = Pilot()
    pilot.resource_name = "cartesius"
    pilot.node_count = 1
    pilot.time = 1 | units.hour
    pilot.queue_name = "short"
    pilot.slots_per_node = 24
    pilot.node_label = "cartesius"
    return pilot


def start_distributed(lgm_node_names):
    lgm_nodes = [new_lgm_node(lgm_node_name) for lgm_node_name in lgm_node_names]
    instance = DistributedAmuse(redirection="file", redirect_file="distributed_amuse.log")
    instance.initialize_code()
    instance.resources.add_resource(new_lgm_gateway())
    for lgm_node in lgm_nodes:
        instance.resources.add_resource(lgm_node)
        
    instance.resources.add_resource(new_cartesius_resource())
    
    instance.pilots.add_pilot(new_local_pilot())
    
    for lgm_node in lgm_nodes:
        instance.pilots.add_pilot(new_gpu_node_pilot(lgm_node))
        instance.pilots.add_pilot(new_cpu_node_pilot(lgm_node))
        
    instance.pilots.add_pilot(new_cartesius_pilot())
    
    print "Pilots:"
    print instance.pilots
    print "Waiting for pilots"
    instance.wait_for_pilots()
    return instance



if __name__ == "__main__":
    logging.basicConfig(filename='example.log', level=logging.WARN)
    logging.getLogger("code").setLevel(logging.DEBUG)
    instance = start_distributed(lgm_node_names=["node06", "node07"])
    try:
        cProfile.run("simulate_young_star_cluster()", "prof")
    finally:
        instance.stop()
        print "Done"
