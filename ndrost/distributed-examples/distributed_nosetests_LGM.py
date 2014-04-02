#!/usr/bin/python
import nose
from amuse.lab import *
from amuse.community.distributed.interface import DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Pilot, Pilots

print "Setting up distributed code"
instance = DistributedAmuse(redirection='none')
#instance.parameters.start_hubs = False

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

pilot = Pilots(2)
pilot.resource_name=['LGMnode00','LGMnode02']
pilot.node_count=1
pilot.time= 2|units.hour
pilot.slots_per_node=1
pilot.node_label=['default', 'default']
print pilot
#~instance.pilots.add_pilot(pilot[0])
instance.pilots.add_pilot(pilot[0])
instance.pilots.add_pilot(pilot[1])
print "Pilots:"
print instance.pilots

print "Waiting for pilots"
instance.wait_for_pilots()

print "Running tests"

nose.run()

print "all tests done, stopping distributed code"

instance.stop()
