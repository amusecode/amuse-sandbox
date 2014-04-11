#!/usr/bin/python
import nose
from amuse.lab import *

from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Pilot, Pilots

#Simple script to run nosetests using the distributed code. Should work for any existing amuse test.
#This example only runs the tests on the local machine.

print "Setting up distributed code"
instance = DistributedAmuse(redirection='none')
instance.initialize_code()

#Add some resources
resource = Resource()
resource.name = "cartesius"
resource.location = "ndrosta@int2-bb.cartesius.surfsara.nl"
resource.amuse_dir = "/home/ndrosta/amuse-svn"
resource.scheduler_type = "slurm"

instance.resources.add_resource(resource)
print "Resources:"
print instance.resources

#Claim nodes on the resources. In this example simply the "local" machine
pilot = Pilot()
pilot.resource_name = "cartesius"
pilot.node_count = 1
pilot.time = 1 | units.hour
pilot.queue_name = "short"
pilot.slots_per_node = 32
pilot.node_label = "cartesius"
instance.pilots.add_pilot(pilot)
print "Pilots:"
print instance.pilots

print "Waiting for pilots"
instance.wait_for_pilots()

print "Running tests"

nose.run()

print "all tests done, stopping distributed code"

instance.stop()
