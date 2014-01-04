"""
 run script for clustergas runs using DistributedAmuse

"""

import sys
import os.path
import numpy
import clustergas_gadget as clustergas
import logging
from amuse.units import units
from amuse.community.fi.interface import Fi
from amuse.community.octgrav.interface import Octgrav
from amuse.community.mi6.interface import MI6
from amuse.community.gadget2.interface import Gadget2
from amuse.community.hermite0.interface import Hermite
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.ph4.interface import ph4
from amuse.community.fastkick.interface import FastKick
from amuse.community.distributed.interface import DistributedAmuse, Resource, Reservation
from SSEplus import SSEplus
import cProfile

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

def new_local_reservation():
    reservation = Reservation()
    reservation.resource_name = "local"
    reservation.node_count = 1
    reservation.time = 2 | units.hour
    reservation.slots_per_node = 10
    reservation.node_label = "local"
    return reservation

def new_cartesius_resource():
    resource = Resource()
    resource.name = "cartesius"
    resource.location = "ndrosta@int2-bb.cartesius.surfsara.nl"
    resource.amuse_dir = "/home/ndrosta/amuse"
    resource.scheduler_type = "slurm"
#  
    return resource


def new_cartesius_reservation():
    reservation = Reservation()
    reservation.resource_name = "cartesius"
    reservation.node_count = 1
    reservation.time = 1 | units.hour
    reservation.queue_name = "short"
    reservation.slots_per_node = 24
    reservation.node_label = "cartesius"
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

    instance.resources.add_resource(new_cartesius_resource())
    
    instance.reservations.add_reservation(new_local_reservation())
    for lgm_node in lgm_nodes:
        instance.reservations.add_reservation(new_gpu_node_reservation(lgm_node))
        instance.reservations.add_reservation(new_cpu_node_reservation(lgm_node))

    instance.reservations.add_reservation(new_cartesius_reservation())
    
    print "Reservations:"
    print instance.reservations
    print "Waiting for reservations"
    instance.wait_for_reservations()
    print "Reservations status after waiting:"
    print instance.reservations
    return instance



if __name__ == "__main__":
    logging.basicConfig(filename='example.log', level=logging.DEBUG)
    #~instance = start_distributed(lgm_node_names=["node12"])#, "node07"])
    instance = start_distributed(lgm_node_names=["node11","node07"])
    numpy.random.seed(123491)
    try:
        kwargs = dict(
            sfeff = 0.3, 
            Nstar = 100,
            Ngas = 100000,
            Rscale = 0.5 | units.parsec,
            runid = "distributed_test",
            feedback_efficiency = 0.01,
            t_end = 0.010 | units.Myr,
            dt_plot = 0.010 | units.Myr,
            #~grav_code = PhiGRAPE,
            grav_code = ph4,
            grav_code_extra = dict(node_label="GPU", mode="gpu", redirection="file", redirect_file="grav_code.log"),
            gas_code = Gadget2,
            gas_code_extra = dict(node_label="cartesius", mode="nogravity", number_of_workers=2, redirection="file", redirect_file="gas_code.log"),
            se_code = SSEplus,
            se_code_extra = dict(redirection="file", redirect_file="se_code.log"),
            #~grav_couple_code = Octgrav,
            #~grav_couple_code_extra = dict(node_label="GPU", redirection="file", redirect_file="grav_couple_code.log")
            grav_couple_code = FastKick,
            grav_couple_code_extra = dict(node_label="GPU", mode="gpu", number_of_workers=1, redirection="file", redirect_file="grav_couple_code.log")
            #~grav_couple_code = BHTree,
            #~grav_couple_code_extra = dict(redirection="file", redirect_file="grav_couple_code.log")
        )
        cProfile.run("clustergas.clustergas(**kwargs)", "prof")
    finally:
        instance.stop()
