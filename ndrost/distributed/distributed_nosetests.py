#!/usr/bin/python
import nose
from amuse.lab import *
from interface import DistributedAmuseInterface
from interface import DistributedAmuse

print "setting up distributed code"
instance = DistributedAmuse(redirection='none')
instance.initialize_code()

resource = Particle()
resource.name='DAS4-VU'
resource.location="niels@fs0.das4.cs.vu.nl"
resource.scheduler_type="sge"
resource.amuse_dir="/home/niels/amuse"
instance.resources.add_particle(resource)
print instance.resources

reservation = Particle()
reservation.resource_name='local'
reservation.node_count=1
reservation.time= 2|units.hour
reservation.slots_per_node=2
reservation.node_label='local'
instance.reservations.add_particle(reservation)

instance.wait_for_reservations()

print "running tests"

nose.run()

print "all tests done, stopping distributed code"

instance.stop()
