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

def new_local_reservation():
    reservation = Reservation()
    reservation.resource_name = "local"
    reservation.node_count = 1
    reservation.time = 2 | units.hour
    reservation.slots_per_node = 10
    reservation.node_label = "local"
    return reservation

def start_distributed():
    instance = DistributedAmuse(redirection="file", redirect_file="distributed_amuse.log")
    instance.initialize_code()
    
    instance.reservations.add_reservation(new_local_reservation())
    
    print "Reservations:"
    print instance.reservations
    print "Waiting for reservations"
    instance.wait_for_reservations()
    return instance


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    instance = start_distributed()
    numpy.random.seed(123491)
    try:
        kwargs = dict(
            sfeff = 0.3, 
            Nstar = 100,
            Ngas = 100000,
            Rscale = 0.5 | units.parsec,
            runid = "distributed_test",
            feedback_efficiency = 0.01,
            t_end = 10.010 | units.Myr,
            dt_plot = 0.010 | units.Myr,
            #~grav_code = PhiGRAPE,
            grav_code = ph4,
            grav_code_extra = dict(node_label="local", redirection="file", redirect_file="grav_code.log"),
            gas_code = Gadget2,
            gas_code_extra = dict(node_label="local", mode="nogravity", number_of_workers=1, redirection="file", redirect_file="gas_code.log"),
            se_code = SSEplus,
            se_code_extra = dict(redirection="file", redirect_file="se_code.log"),
            #~grav_couple_code = Octgrav,
            #~grav_couple_code_extra = dict(node_label="GPU", redirection="file", redirect_file="grav_couple_code.log")
            grav_couple_code = BHTree,
            grav_couple_code_extra = dict(node_label="local", redirection="file", redirect_file="grav_couple_code.log")
        )
        cProfile.run("clustergas.clustergas(**kwargs)", "prof")
    finally:
        instance.stop()
