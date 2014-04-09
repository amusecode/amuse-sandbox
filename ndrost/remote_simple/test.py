from amuse.community.distributed.interface import DistributedAmuse
from amuse.community.distributed.distributed_datamodel import Pilot
from simple import simpleInterface,Simple

from amuse.lab import *

distinstance = DistributedAmuse(redirection='none')
distinstance.parameters.debug = True

print "Resources:"
print distinstance.resources

pilot = Pilot()
pilot.resource_name='local'
pilot.node_count=1
pilot.time= 2|units.hour
pilot.slots_per_node=2
pilot.label='local'
distinstance.pilots.add_pilot(pilot)
print "Pilots:"
print distinstance.pilots

print "Waiting for pilots"
distinstance.wait_for_pilots()
distinstance.set_as_default_for_all_workers()

s=Simple(redirection="none", dynamic_python_code=True)

print s.timestwo(5.)


s.stop()
