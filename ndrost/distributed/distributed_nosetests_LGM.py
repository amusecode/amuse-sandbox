#!/usr/bin/python
import nose
from amuse.lab import *
from interface import DistributedAmuseInterface
from interface import DistributedAmuse
from interface import Resource, Resources, Reservation, Reservations

print "Setting up distributed code"
instance = DistributedAmuse(redirection='none')
instance.initialize_code()

resource = Resource()
resource.name='LGM'
#~resource.location="vriesn@node00"
#~resource.gateway="vriesn@fs.lgm.liacs.nl"
resource.location="fs.lgm.liacs.nl"
resource.amuse_dir="/home/niels/amuse"
instance.resources.add_resource(resource)

resource = Resource()
resource.name="LGMnode00"
resource.location="node00"
resource.gateway="LGM"
resource.amuse_dir="/home/niels/amuse"
instance.resources.add_resource(resource)
print "Resources:"
print instance.resources

reservation = Reservations(2)
reservation.resource_name=['local', 'LGMnode00']
reservation.node_count=1
reservation.time= 2|units.hour
reservation.slots_per_node=10
reservation.node_label=['local', 'LGMnode00']
print reservation
#~instance.reservations.add_reservation(reservation[0])
instance.reservations.add_reservation(reservation[1])
print "Reservations:"
print instance.reservations

print "Waiting for reservations"
instance.wait_for_reservations()

print "Running tests"

nose.run()

print "all tests done, stopping distributed code"

instance.stop()
