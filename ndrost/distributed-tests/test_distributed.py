
from amuse.community.distributed.interface import DistributedAmuse, Pilot, Resource, Pilots, Resources

from amuse.test.amusetest import TestCase
from amuse.support import exceptions

from amuse.units import units

import subprocess
import os
import time
import unittest
import logging

logger = logging.getLogger(__name__)

class TestDistributedLocal(TestCase):
    
    def initialize_distributed_code(self):
        logger.info("initializing distributed code")
        distinstance = self.create_distributed_code()
        
        self.create_resources(distinstance)
        
        self.create_pilot(distinstance, True)
        
        return distinstance
    
    def create_distributed_code(self):
        logger.info("Setting up distributed code")
        distinstance = DistributedAmuse()#redirection='none')
        distinstance.parameters.debug = False
        
        return distinstance
    
    def create_resources(self, distinstance):
        pass
    
    def create_pilot(self, distinstance, wait=True):
        logger.info("starting local pilot")
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
        "Test if default local resource is created"

        distinstance = self.create_distributed_code()
        
        distinstance.commit_parameters()
        
        local_resource = distinstance.resources[0]
        
        self.assertEqual(local_resource.name, "local")
        
        distinstance.stop()
        
    def test2b(self):
        "Test creating default resource"

        distinstance = self.create_distributed_code()
        
        self.create_resources(distinstance)
        
        distinstance.stop()

    def test2c(self):
        "Test creating additional (local) resource"

        distinstance = self.create_distributed_code()
        
        resource = Resource()
        resource.name = "also-local"
        resource.location = "local"
        resource.amuse_dir = self.get_amuse_root_dir()
        resource.scheduler_type = "local"
        
        distinstance.resources.add_resource(resource)
        
        distinstance.stop()
        
        
    def test2d(self):
        "Test creating and removing additional (local) resource"

        distinstance = self.create_distributed_code()
        
        resource = Resource()
        resource.name = "also-local"
        resource.location = "local"
        resource.amuse_dir = self.get_amuse_root_dir()
        resource.scheduler_type = "local"

        distinstance.resources.add_resource(resource)
        
        logger.info(distinstance.resources)
        
        self.assertEquals(len(distinstance.resources), 2)
        
        distinstance.resources.remove_resource(resource)
        
        self.assertEquals(len(distinstance.resources), 1)
        
        distinstance.stop()
        
        
    def test2e(self):
        "Test creating multiple (local) resources"

        distinstance = self.create_distributed_code()

        resources = Resources()
        for i in range(1, 11):       
            resource = Resource()
            resource.name = "also-local-" + str(i)
            resource.location = "local"
            resource.amuse_dir = self.get_amuse_root_dir()
            resource.scheduler_type = "local"

            resources.add_resource(resource)
        
        distinstance.resources.add_resources(resources)
        
        logger.info(distinstance.resources)
        
        #one default local resources, and 10 additional resources created above
        self.assertEquals(len(distinstance.resources), 11)
        
        self.assertEquals(distinstance.resources[0].name, "local")
        
        for i in range(1, 11):
            self.assertEquals(distinstance.resources[i].name, "also-local-" + str(i))
        
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
        
        distinstance = self.initialize_distributed_code()

        self.assertEquals(len(distinstance.pilots), 1)
        
        distinstance.pilots.remove_pilot(distinstance.pilots[0])
        
        self.assertEquals(len(distinstance.pilots), 0)
        
        distinstance.stop
        
class TestDistributedCartesius(TestDistributedLocal):
    
    def create_resources(self, distinstance):
        logger.info("creating cartesius resource")
        resource = Resource()
        resource.name = "cartesius"
        resource.location = "ndrosta@int2-bb.cartesius.surfsara.nl"
        resource.amuse_dir = "/home/ndrosta/amuse"
        resource.scheduler_type = "slurm"
        distinstance.resources.add_resource(resource)
        logger.info("creating cartesius resource done")

    
    def create_pilot(self, distinstance, wait=True):
        logger.info("starting cartesius pilot")
        pilot = Pilot()
        pilot.resource_name = "cartesius"
        pilot.node_count = 1
        pilot.time = 1 | units.hour
        pilot.queue_name = "short"
        pilot.slots_per_node = 10
        pilot.label = "cartesius"
        distinstance.pilots.add_pilot(pilot)

        if wait:
            logger.info("Waiting for cartesius pilot")
            distinstance.wait_for_pilots()
            logger.info("Done waiting for cartesius pilot")
            
        
#     def test3c(self):
#         self.skip("this test specific for local machine")
        
        
