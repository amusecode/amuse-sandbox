#!/usr/bin/python
import sys
from amuse.lab import *
from amuse.support.exceptions import AmuseException

from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Reservation, Reservations


print "Setting up distributed code"
instance = DistributedAmuse(redirection='none')
#instance = DistributedAmuse()
instance.initialize_code()

try:
    resource = Resource()
    resource.name='DAS4-VU'
    resource.location="niels@fs0.das4.cs.vu.nl"
    resource.scheduler_type="sge"
    resource.amuse_dir="/home/niels/amuse"
    instance.resources.add_resource(resource)
    print "Resources:"
    print instance.resources
    
    resource = Resource()
    resource.name='DAS4-Leiden'
    resource.location="niels@fs1.das4.liacs.nl"
    resource.scheduler_type="sge"
    resource.amuse_dir="/home/niels/amuse"
    instance.resources.add_resource(resource)
    print "Resources:"
    print instance.resources
    
    reservation = Reservation()
    reservation.resource_name='DAS4-Leiden'
    reservation.node_count=1
    reservation.time= 2|units.hour
    reservation.slots_per_node=22
#    reservation.queue_name="all.q"
    reservation.node_label='normal'
#    reservation.options="resources=gpu=GTX480"
    instance.reservations.add_reservation(reservation)
    
    
    print "Reservations:"
    print instance.reservations

    print "Waiting for reservations"
    instance.wait_for_reservations()

except AmuseException:
    print "Error while obtaining resources:", instance.get_current_error()
    raise

print "Running script"

script = sys.argv[1]

sys.argv = sys.argv[1:]

execfile(script)

print "all tests done, stopping distributed code"

instance.stop()
