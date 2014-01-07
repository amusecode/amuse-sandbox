"""
 run script for clustergas runs using DistributedAmuse

"""

import logging
import cProfile

from amuse.community.distributed.interface import DistributedAmuse, Resource, Reservation

from young_star_cluster import simulate_young_star_cluster

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

def new_local_reservation():
    reservation = Reservation()
    reservation.resource_name = "local"
    reservation.node_count = 1
    reservation.time = 2 | units.hour
    reservation.slots_per_node = 10
    reservation.node_label = "local"
    return reservation

def new_gpu_node_reservation(resource):
    reservation = Reservation()
    reservation.resource_name = resource.name
    reservation.node_count = 1
    reservation.time = 2 | units.hour
    reservation.slots_per_node = 2
    reservation.node_label = "GPU"
    return reservation

def new_cpu_node_reservation(resource):
    reservation = Reservation()
    reservation.resource_name = resource.name
    reservation.node_count = 1
    reservation.time = 2 | units.hour
    reservation.slots_per_node = 10
    reservation.node_label = "CPU"
    return reservation


def start_distributed(lgm_node_names):
    lgm_nodes = [new_lgm_node(lgm_node_name) for lgm_node_name in lgm_node_names]
    instance = DistributedAmuse(redirection="file", redirect_file="distributed_amuse.log")
    instance.initialize_code()
    instance.resources.add_resource(new_lgm_gateway())
    for lgm_node in lgm_nodes:
        instance.resources.add_resource(lgm_node)
    
    instance.reservations.add_reservation(new_local_reservation())
    for lgm_node in lgm_nodes:
        instance.reservations.add_reservation(new_gpu_node_reservation(lgm_node))
        instance.reservations.add_reservation(new_cpu_node_reservation(lgm_node))
    
    print "Reservations:"
    print instance.reservations
    print "Waiting for reservations"
    instance.wait_for_reservations()
    return instance



if __name__ == "__main__":
    logging.basicConfig(filename='example.log', level=logging.WARN)
    logging.getLogger("code").setLevel(logging.DEBUG)
    #~instance = start_distributed(lgm_node_names=["node12"])#, "node07"])
#~    instance = start_distributed(lgm_node_names=["node12", "node07", "node03"])
    try:
        cProfile.run("simulate_young_star_cluster()", "prof")
    finally:
#~        instance.stop()
        print "Done"
