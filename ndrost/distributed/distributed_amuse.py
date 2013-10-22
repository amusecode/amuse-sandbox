#!/usr/bin/python
import sys
from amuse.lab import *
from interface import DistributedAmuseInterface
from interface import DistributedAmuse
from interface import Resource, Resources, Reservation, Reservations

print "Setting up distributed code"
instance = DistributedAmuse(redirection='none')
instance.initialize_code()

resource = Resource()
#resource.name='DAS4-VU'
#resource.location="niels@fs0.das4.cs.vu.nl"
#resource.scheduler_type="sge"
#resource.amuse_dir="/home/niels/amuse"
#instance.resources.add_resource(resource)
print "Resources:"
print instance.resources

reservation = Reservation()
reservation.resource_name='local'
reservation.node_count=1
reservation.time= 2|units.hour
reservation.slots_per_node=22
reservation.node_label='local'
instance.reservations.add_reservation(reservation)

print "Reservations:"
print instance.reservations

print "Waiting for reservations"
instance.wait_for_reservations()

print "Running script"

script = sys.argv[1]

sys.argv = sys.argv[1:]

execfile(script)

print "all tests done, stopping distributed code"


