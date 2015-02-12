import sys
from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Pilots, Pilot
from amuse.units import units

def new_lgm_gateway():
    resource = Resource()
    resource.name = "LGM"
    resource.location = "fs.lgm.liacs.nl"
    resource.amuse_dir = "/scratch/jilkova/amuse/amuse-svn"
    return resource

def new_lgm_node(node_name):
    resource = Resource()
    resource.name = "LGM" + node_name
    resource.location = node_name
    resource.gateway = "LGM"
    resource.amuse_dir = "/scratch/jilkova/amuse/amuse-svn"
    return resource

def new_local_pilot():
    pilot = Pilot()
    pilot.resource_name = "local"
    pilot.node_count = 1
    #pilot.time = 2 | units.hour
    pilot.time = 666 | units.day
    pilot.slots_per_node = 12
    pilot.node_label = "local"
    return pilot

def new_cpu_node_pilot(resource):
    pilot = Pilot()
    pilot.resource_name = resource.name
    pilot.node_count = 1
    #pilot.time = 2 | units.hour
    pilot.time = 666 | units.day
    pilot.slots_per_node = 12
    pilot.node_label = "CPU"
    return pilot
  
def start_distributed(lgm_node_names):
    lgm_nodes = [new_lgm_node(lgm_node_name) for lgm_node_name in lgm_node_names]
#   instance = DistributedAmuse()
    instance = DistributedAmuse(redirection="file", redirect_file="distributed_amuse.log")
#   instance = DistributedAmuse(redirection="none")
    instance.initialize_code()
    instance.use_for_all_workers(True)
    instance.resources.add_resource(new_lgm_gateway())
    for lgm_node in lgm_nodes:
        instance.resources.add_resource(lgm_node)
    
    instance.pilots.add_pilot(new_local_pilot())
    for lgm_node in lgm_nodes:
        instance.pilots.add_pilot(new_cpu_node_pilot(lgm_node))
    
    print " ** distributed amuse ** "
    print " ** pilots:"
    print instance.pilots
    print "\t waiting for pilots"
    instance.wait_for_pilots()
    return instance

if __name__ == "__main__":
    #logging.basicConfig(filename='example.log', level=logging.WARN)
    #logging.getLogger("code").setLevel(logging.DEBUG)
    
    instance = start_distributed(lgm_node_names=["node17", "node16"])
    try:
        mpi_test_2()
    finally:
        instance.stop()
        print " ** done "

