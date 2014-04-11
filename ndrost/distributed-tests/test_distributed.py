
from amuse.community.distributed.interface import DistributedAmuse, Pilot, Resource, Pilots

from amuse.test.amusetest import TestCase
from amuse.support import exceptions

from amuse.units import units

import subprocess
import os
import time
import unittest


class TestDistributedLocal(TestCase):
    
    def initialize_distributed_code(self):
        "Default initialization of distributed code for all tests. Creates a single pilot"
        distinstance = self.create_distributed_code()
        
        self.create_resources(distinstance)
        
        self.create_pilot(distinstance, True)
        
        return distinstance
    
    def create_distributed_code(self):
        print "Setting up distributed code"
        distinstance = DistributedAmuse(redirection='none')
        distinstance.parameters.debug = True
        
        return distinstance
    
    def create_resources(self, distinstance):
        pass
    
    def create_pilot(self, distinstance, wait=True):
        print "starting local pilot"
        pilot = Pilot()
        pilot.resource_name='local'
        pilot.node_count=1
        pilot.time= 2|units.hour
        pilot.slots_per_node=10
        pilot.label='local'
        distinstance.pilots.add_pilot(pilot)
        
        if wait:
            distinstance.wait_for_pilots()
    
    #CODE
    
    def test1a(self):
        "Basic test of init of distributed code"

        distinstance = self.create_distributed_code()
        
        distinstance.stop()

    #RESOURCES
        
    def test2a(self):
        "Test creating default resource"

        distinstance = self.create_distributed_code()
        
        self.create_resources(distinstance)
        
        distinstance.stop()

    def test2b(self):
        "Test creating additional (local) resource"

        distinstance = self.create_distributed_code()
        
        resource = Resource()
        resource.name = "also-local"
        resource.location = "local"
        resource.amuse_dir = self.get_amuse_root_dir()
        resource.scheduler_type = "local"
        
        distinstance.resources.add_resource(resource)
        
        distinstance.stop()
        
        
    def test2c(self):
        "Test creating and removing additional (local) resource"

        distinstance = self.create_distributed_code()
        
        resource = Resource()
        resource.name = "also-local"
        resource.location = "local"
        resource.amuse_dir = self.get_amuse_root_dir()
        resource.scheduler_type = "local"
        
        distinstance.resources.add_resource(resource)
        
        distinstance.resources.remove(resource)
        
        distinstance.stop()

        
    
    

    #PILOTS
    
    def test3a(self):
        "Test creating resources and starting a pilot"

        distinstance = self.initialize_distributed_code()
        
        distinstance.stop()
        
    def test3b(self):
        "Test starting multiple pilots one by one"

        distinstance = self.create_distributed_code()
        self.create_resources(distinstance)

        for i in range(0, 10):
            self.create_pilot(distinstance, wait=False)
        
        distinstance.wait_for_pilots()
        
        distinstance.stop()
        
    def test3c(self):
        "Test starting multiple (local) pilots at once"

        distinstance = self.create_distributed_code()

        pilots = Pilots()

        for i in range(0, 10):        
            pilot = Pilot()
            pilot.resource_name='local'
            pilot.node_count=1
            pilot.time= 2|units.hour
            pilot.slots_per_node=10
            pilot.label='local'
            pilots.add_pilot(pilot)
        
        distinstance.pilots.add_pilots(pilots)
        
        distinstance.wait_for_pilots()
        
        distinstance.stop()
        
    def test3d(self):
        "Cancel and remove a pilot"
        
        self.skip("this should work, but doesn't")
        
        distinstance = self.initialize_distributed_code()
        
        print distinstance.pilots
        
        #distinstance.pilots.remove_particle(distinstance.pilots[0])
        distinstance.pilots.delete_pilot(distinstance.pilots[0])
        
        print distinstance.pilots
        
        distinstance.stop
        


class TestDistributedCartesius(TestDistributedLocal):
    
    def create_resources(self, distinstance):
        resource = Resource()
        resource.name = "cartesius"
        resource.location = "ndrosta@int2-bb.cartesius.surfsara.nl"
        resource.amuse_dir = "/home/ndrosta/amuse"
        resource.scheduler_type = "slurm"
        distinstance.resources.add_resource(resource)
    
    def create_pilot(self, distinstance, wait=True):
        print "starting cartesius pilot"
        pilot = Pilot()
        pilot.resource_name = "cartesius"
        pilot.node_count = 1
        pilot.time = 1 | units.hour
        pilot.queue_name = "short"
        pilot.slots_per_node = 10
        pilot.label = "cartesius"
        distinstance.pilots.add_pilot(pilot)

        if wait:
            print "Waiting for pilots"
            distinstance.wait_for_pilots()
        
    def test4(self):
        self.skip("this test specific for local machine")
        
        