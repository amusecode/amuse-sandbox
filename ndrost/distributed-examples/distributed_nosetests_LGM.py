#!/usr/bin/python
import nose
from amuse.lab import *
from amuse.community.distributed.interface import DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Reservation, Reservations

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

resource2 = Resource()
resource2.name="LGMnode02"
resource2.location="node02"
resource2.gateway="LGM"
resource2.amuse_dir="/home/niels/amuse"
instance.resources.add_resource(resource2)
print "Resources:"
print instance.resources

print instance.resources

reservation = Reservations(2)
reservation.resource_name=['LGMnode00','LGMnode02']
reservation.node_count=1
reservation.time= 2|units.hour
reservation.slots_per_node=1
reservation.node_label=['default', 'default']
print reservation
#~instance.reservations.add_reservation(reservation[0])
instance.reservations.add_reservation(reservation[0])
instance.reservations.add_reservation(reservation[1])
print "Reservations:"
print instance.reservations

print "Waiting for reservations"
instance.wait_for_reservations()

print "Running tests"

nose.run()

print "all tests done, stopping distributed code"

instance.stop()
