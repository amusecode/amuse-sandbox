
from amuse.community.distributed.interface import DistributedAmuse, Pilot, Resource, Pilots, Resources

from amuse.test.amusetest import TestCase
from amuse.support import exceptions, options

from amuse.units import units

from amuse.rfi.channel import DistributedChannel

from amuse.community.bhtree.interface import BHTree

import subprocess
import os
import time
import unittest
import logging
import socket
import hashlib

logger = logging.getLogger(__name__)

from amuse.community import *

from amuse.rfi.core import PythonCodeInterface

from amuse.support.interface import InCodeComponentImplementation

class UtilityImplementation(object):

  @staticmethod
  def caclulate_hash(filename):
      hasher = hashlib.sha256()

      file = open(filename, 'rb')

      filecontent = file.read()

      file.close()

      hasher.update(filecontent)

      return hasher.hexdigest()

  def timestwo(self, xin, xout):
    print "timestwo is running", os.getcwd()

    xout.value = 2 * xin
    return 0

  def hash(self, filename, hash):
    print "calculating hash of", filename

    hash.value = UtilityImplementation.caclulate_hash(filename)

    return 0

  def hostname(self, hostname):
    print "returning hostname", socket.gethostname()

    hostname.value = socket.gethostname()
    return 0

  def cwd(self, cwd):
    print "returning cwd ", os.getcwd()

    cwd.value = os.getcwd()
    return 0




class UtilityInterface(PythonCodeInterface):

    def __init__(self, **options):
        PythonCodeInterface.__init__(self, UtilityImplementation, **options)

    @legacy_function
    def timestwo():
        function = LegacyFunctionSpecification()
        function.addParameter('xin', dtype='float64', direction=function.IN, unit=units.none)
        function.addParameter('xout', dtype='float64', direction=function.OUT, unit=units.km)
        function.result_type = 'int32'
        return function

    @legacy_function
    def hash():
        function = LegacyFunctionSpecification()
        function.addParameter('filename', dtype='string', direction=function.IN)
        function.addParameter('hash', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def hostname():
        function = LegacyFunctionSpecification()
        function.addParameter('hostname', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def cwd():
        function = LegacyFunctionSpecification()
        function.addParameter('cwd', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function



class Utility(InCodeComponentImplementation):
    def __init__(self, convert_nbody=None, **options):
        interface = UtilityInterface(**options)

        InCodeComponentImplementation.__init__(self, interface, **options)



class TestDistributedLocal(TestCase):

    def initialize_distributed_code(self):
        logger.info("initializing distributed code")
        distinstance = self.create_distributed_code()

        self.create_resources(distinstance)

        self.create_pilot(distinstance, True)

        return distinstance

    def create_distributed_code(self):
        logger.info("Setting up distributed code")
        distinstance = DistributedAmuse(redirection='none')
        distinstance.parameters.debug = True

        return distinstance

    def create_resources(self, distinstance):
        pass

    def create_pilot(self, distinstance, wait=True):
        logger.info("starting local pilot")
        pilot = Pilot()
        pilot.resource_name = 'local'
        pilot.node_count = 1
        pilot.time = 2 | units.hour
        pilot.slots_per_node = 10
        pilot.label = 'local'
        distinstance.pilots.add_pilot(pilot)

        if wait:
            distinstance.wait_for_pilots()

    # CODE

    def test1a(self):
        "Basic test of init of distributed code"

        distinstance = self.create_distributed_code()

        distinstance.stop()

    # RESOURCES

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

        # one default local resources, and 10 additional resources created above
        self.assertEquals(len(distinstance.resources), 11)

        self.assertEquals(distinstance.resources[0].name, "local")

        for i in range(1, 11):
            self.assertEquals(distinstance.resources[i].name, "also-local-" + str(i))

        distinstance.stop()

    # PILOTS

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
            pilot.resource_name = 'local'
            pilot.node_count = 1
            pilot.time = 2 | units.hour
            pilot.slots_per_node = 10
            pilot.label = 'local'
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

        distinstance.stop()


    # WORKER JOBS

    def test4a(self):
        "Test setting code as default for all workers"

        distinstance = self.initialize_distributed_code()

        distinstance.use_for_all_workers()

        # check if distributed channel is now set as default
        self.assertEqual(DistributedChannel.default_distributed_instance, distinstance)
        self.assertTrue(options.GlobalOptions.instance().overriden_options.has_key("channel_type"))
        self.assertEquals(options.GlobalOptions.instance().overriden_options["channel_type"], "distributed")

        distinstance.stop()


    def test4b(self):
        "Test creating distributed worker with explicit distributed code"

        distinstance = self.initialize_distributed_code()

        self.assertEqual(len(distinstance.workers), 0)

        bhtree = BHTree(channel_type="distributed", distributed_instance=distinstance)

        bhtree.stop()

        # check if the distributed instance has a worker (and thus the setting had an effect)
        self.assertEqual(len(distinstance.workers), 1)

        distinstance.stop()


    def test4c(self):
        "Test creating distributed worker with default distributed code"

        distinstance = self.initialize_distributed_code()

        distinstance.use_for_distributed_workers()

        self.assertEqual(len(distinstance.workers), 0)

        bhtree = BHTree(channel_type="distributed")

        bhtree.stop()

        # check if the distributed instance has a worker (and thus the setting had an effect)
        self.assertEqual(len(distinstance.workers), 1)

        distinstance.stop()

    def test4d(self):
        "Test creating distributed worker with distributed code implcitly"

        distinstance = self.initialize_distributed_code()

        distinstance.use_for_all_workers()

        self.assertEqual(len(distinstance.workers), 0)

        bhtree = BHTree()

        bhtree.stop()

        # check if the distributed instance has a worker (and thus the setting had an effect)
        self.assertEqual(len(distinstance.workers), 1)

        distinstance.stop()

    def test4e(self):
        "Test if creating a distributed code with any but the sockets channel leads to an error"

        distinstance = DistributedAmuse()
        distinstance.commit_parameters()
        distinstance.use_for_all_workers()

        try:
            distinstance2 = DistributedAmuse()
        except Exception:
            # exception we are hoping for
            return
        finally:
            distinstance.stop()


        raise Exception("Test failed, as distributed channel did not produce an error")


    def test4f(self):
        "Test if creating two distributed codes and using one of them works"

        distinstance1 = self.initialize_distributed_code()

        distinstance2 = self.initialize_distributed_code()

        distinstance1.use_for_all_workers()

        distinstance2.use_for_all_workers()

        self.assertEqual(len(distinstance2.workers), 0)

        bhtree = BHTree()

        bhtree.stop()

        # check if the correct distributed instance has a worker (and thus the setting had an effect)
        self.assertEqual(len(distinstance2.workers), 1)

        distinstance1.stop()
        distinstance2.stop()

    # DYNAMIC PYTHON WORKER JOBS

    def test5a(self):
        "Test starting a dynamic python worker"

        distinstance = self.initialize_distributed_code()

        distinstance.use_for_all_workers()

        u = Utility(redirection="none", dynamic_python_code=True)

        print u.timestwo(5.)

        u.stop()

        distinstance.stop()

    def test5b(self):
        "Test if files get copied to the remote site properly"

        distinstance = self.initialize_distributed_code()

        distinstance.use_for_all_workers()

        u = Utility(redirection="none", dynamic_python_code=True)

        print u.timestwo(5.)

        for file in ['test_distributed.py', 'bigfile.py']:
            remote_hash = u.hash(file)[0]
            local_hash = UtilityImplementation.caclulate_hash(file)
            self.assertEqual(remote_hash, local_hash)

        print u.hostname()
        print u.cwd()

        u.stop()

        distinstance.stop()

    def test5c(self):
        "Test if setting the tmp dir of a resource has an effect"
        
        test_tmp_dir = os.path.expanduser("~/distributed-test-5c")

        distinstance = self.create_distributed_code()

        resource = Resource()
        resource.name = "also-local"
        resource.location = "local"
        resource.amuse_dir = self.get_amuse_root_dir()
        resource.tmp_dir = test_tmp_dir
        resource.scheduler_type = "local"
        distinstance.resources.add_resource(resource)

        pilot = Pilot()
        pilot.resource_name = 'also-local'
        pilot.node_count = 1
        pilot.time = 2 | units.hour
        pilot.slots_per_node = 10
        pilot.label = 'also-local'
        distinstance.pilots.add_pilot(pilot)

        distinstance.wait_for_pilots()

        distinstance.use_for_all_workers()

        u = Utility(redirection="none", dynamic_python_code=True)

        cwd = u.cwd()

        u.stop()

        distinstance.stop()

        print "cwd while running test was", cwd
        
        self.assertTrue(cwd.startswith(test_tmp_dir))

    # SCRIPT JOBS

    # FUNCTION JOBS


class TestDistributedCartesius(TestDistributedLocal):

    def create_resources(self, distinstance):
        logger.info("creating cartesius resource")
        resource = Resource()
        resource.name = "cartesius"
        resource.location = "ndrosta@int2-bb.cartesius.surfsara.nl"
        resource.amuse_dir = "/home/ndrosta/amuse"
        resource.tmp_dir = "/home/ndrosta/amuse-tmp"
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

    def test5c(self):
        "Test if setting the tmp dir for a cartesius node has an effect"
        
        distinstance = self.initialize_distributed_code()

        distinstance.use_for_all_workers()

        u = Utility(redirection="none", dynamic_python_code=True)

        cwd = u.cwd()

        u.stop()

        distinstance.stop()

        print "cwd while running test was", cwd
        
        self.assertTrue(cwd.startswith('/home/ndrosta/amuse-tmp'))

    def test_cartesius_5d(self):
        "Test if starting a dynamic python worker gets us to a cartesius node"

        distinstance = self.initialize_distributed_code()

        distinstance.use_for_all_workers()

        u = Utility(redirection="none", dynamic_python_code=True)

        hostname = u.hostname()

        u.stop()

        distinstance.stop()

        print "test ran on", hostname

        self.assertTrue(hostname.startswith("tcn"))

