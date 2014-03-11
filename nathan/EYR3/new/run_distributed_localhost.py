"""
 run script for clustergas runs using DistributedAmuse

"""

import logging
import cProfile

from amuse.community.distributed.interface import DistributedAmuse, Resource, Pilot
from amuse.lab import *

from young_star_cluster import simulate_young_star_cluster

def new_local_pilot():
    pilot = Pilot()
    pilot.resource_name = "local"
    pilot.node_count = 1
    pilot.time = 2 | units.hour
    pilot.slots_per_node = 10
    pilot.node_label = "local"
    return pilot

def start_distributed():
    instance = DistributedAmuse(redirection="file", redirect_file="distributed_amuse.log")
    instance.initialize_code()
    
    instance.pilots.add_pilot(new_local_pilot())
    
    print "Pilots:"
    print instance.pilots
    print "Waiting for pilots"
    instance.wait_for_pilots()
    return instance

if __name__ == "__main__":
    logging.basicConfig(filename='run_distributed_localhost.log', level=logging.WARN)
    logging.getLogger("code").setLevel(logging.DEBUG)
    instance = start_distributed()
    simulate_young_star_cluster()
