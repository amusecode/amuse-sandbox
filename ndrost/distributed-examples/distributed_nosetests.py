#!/usr/bin/python
import nose



from amuse.lab import *

from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Pilot, Pilots

#Simple script to run nosetests using the distributed code. Should work for any existing amuse test.
#This example only runs the tests on the local machine.


def new_lgm_gateway():
    resource = Resource()
    resource.name = "LGM"
    resource.location = "fs.lgm.liacs.nl"
    resource.amuse_dir = "amuse"
    return resource

def new_lgm_node(node_name):
    resource = Resource()
    resource.name = "LGM" + node_name
    resource.location = node_name
    resource.gateway = "LGM"
    resource.amuse_dir = "amuse"
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
    resource.amuse_dir = "amuse-svn"
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
    instance.parameters.debug = True
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

    instance = start_distributed(lgm_node_names=["node07", "node08"])

    print "Running tests"

    nose.run()

    print "all tests done, stopping distributed code"

    instance.stop()
