#!/usr/bin/python
import sys
from amuse.units import units
from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Pilot, Pilots
import atexit

class local_only(object):
  def __init__(self):
  
    print "Setting up distributed code"
    instance = DistributedAmuse()
    instance.initialize_code()

#    instance.parameters.debug = True
#    instance.parameters.webinterface_port = 4556
    instance.commit_parameters()
  
    print "url:", instance.get_webinterface_url()
  
    print "Resources:"
    print instance.resources
      
    #Claim nodes on the resources. In this example simply the "local" machine
    pilot = Pilot()
    pilot.resource_name='local' # label of the resource to be used
    pilot.node_count=1 # desired number of nodes
    pilot.time= 99|units.hour # wallclock that resource remains available (mainly for job queues)
    pilot.slots_per_node=99 # slots is accounting measure for a job 
    pilot.label='local' # label for subgroups of the resource 
    
    instance.pilots.add_pilot(pilot)

    print "Reservations:"
    print instance.pilots
    
    print "Waiting for reservations"
    instance.wait_for_pilots()
    self.instance=instance
  def __enter__(self):
    self.instance.set_as_default_for_all_workers()
    return self.instance
  def __exit__(self,exc_type, exc_val, exc_tb):
    self.instance.stop()
    return False  

class local_and_paddegat(object):
  def __init__(self,node_label=None):
  
    print "Setting up distributed code"
    instance = DistributedAmuse()
    instance.initialize_code()

    instance.parameters.debug = False
#    instance.parameters.webinterface_port = 4556
    instance.commit_parameters()
  
    print "url:", instance.get_webinterface_url()
  
  
    #Add some resources
    resource = Resource()
    resource.name='paddegat'
    resource.location="pelupes@paddegat.strw.leidenuniv.nl"
    #resource.scheduler_type="sge"
    resource.amuse_dir="/disks/paddegat2/pelupes/amuse/amuse-strw"
    instance.resources.add_resource(resource)
  
    print "Resources:"
    print instance.resources
      
    #Claim nodes on the resources. In this example simply the "local" machine
    pilot = Pilot()
    pilot.resource_name='local' # label of the resource to be used
    pilot.node_count=1 # desired number of nodes
    pilot.time= 99|units.hour # wallclock that resource remains available (mainly for job queues)
    pilot.slots_per_node=99 # slots is accounting measure for a job 
    pilot.label='local' # label for subgroups of the resource 
    
    instance.pilots.add_pilot(pilot)

    #Claim nodes on the resources. In this example simply the "local" machine
    pilot = Pilot()
    pilot.resource_name='paddegat' # label of the resource to be used
    pilot.node_count=1 # desired number of nodes
    pilot.time= 99|units.hour # wallclock that resource remains available (mainly for job queues)
    pilot.slots_per_node=99 # slots is accounting measure for a job 
    pilot.label='paddegat' # label for subgroups of the resource 
    
    instance.pilots.add_pilot(pilot)

    print "Pilots:"
    print instance.pilots
    
    print "Waiting for pilots"
    instance.wait_for_pilots()
    print "done"
    self.instance=instance
  def __enter__(self):
    self.instance.set_as_default_for_all_workers()
    return self.instance
  def __exit__(self,exc_type, exc_val, exc_tb):
    return False

