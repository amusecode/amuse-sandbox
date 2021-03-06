from amuse.community import *
from amuse.test.amusetest import TestWithMPI
from amuse.lab import *

from amuse.community.distributed.interface import DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Reservation, Reservations

from amuse.community.hermite0.interface import Hermite

class DistributedTests(TestWithMPI):
    
    def start_nodes(self):
	
        instance = DistributedAmuse(redirection='none')
        instance.initialize_code()
        
        resource = Particle()
        resource.name='DAS4-VU'
        resource.location="niels@fs0.das4.cs.vu.nl"
        resource.scheduler_type="sge"
        resource.amuse_dir="/home/niels/amuse"
        instance.resources.add_particle(resource)
        print instance.resources

#        instance.new_resource(name='DAS4-VU',
#                              location="niels@fs0.das4.cs.vu.nl",
#                              scheduler_type="sge", 
#                              amuse_dir="/home/niels/amuse",
#                              )
        
#        instance.new_resource(name='DAS4-Leiden',
#                              hostname="fs1.das4.liacs.nl",
#                              username="niels",
#                              scheduler_type="sge", 
#                              amuse_dir="/home/niels/amuse",
#                              )

 #        instance.new_resource(name='LGM',
#                              hostname="fs.lgm.liacs.nl", 
#                              amuse_dir='/home/niels/amuse-svn',
#                              scheduler_type="local", 
#                              username="niels")
        
#        instance.new_resource(name='LGM-4', 
#                              hostname="node004", 
#                              gateways="fs.lgm.liacs.nl",
#                              scheduler_type="local", 
#                              amuse_dir='/var/local/amuse',
#                              username="niels")
    
#        instance.new_reservation(resource_name='DAS4-VU', node_count=5, time= 2|units.hour, node_label='VU')
#        instance.new_reservation(resource_name='DAS4-Leiden', node_count=5, time= 2|units.hour, node_label='Leiden')
#        instance.new_reservation(resource_name='local', node_count=1, time= 2|units.hour, slots=2, node_label='local')
        
        reservation = Particle()
        reservation.resource_name='local'
        reservation.node_count=1
        reservation.time= 2|units.hour
        reservation.slots=2
        reservation.node_label='local'
        instance.reservations.add_particle(reservation)
#        instance.new_reservation(resource_name='LGM-4', node_count=1, time=2|units.hour, node_label='LGM')
    
        instance.wait_for_reservations()

        return instance

    def test0(self):
        print "starting"
        instance = self.start_nodes()

        print "taking a nap"
        import time
        time.sleep(60)

        print "stopping instance"
        instance.stop()

    def test1(self):
        print "starting"
        instance = self.start_nodes()

        print "starting codes"
        gravity = Hermite(number_of_workers = 1, redirection='none')
        gravity2 = Hermite(number_of_workers = 1, redirection='none')
        
        #gadget2 = Gadget(nr_of_workers=4, nr_of_nodes=2, node_label='VU')

        # some interesting simulation using these workers


        print "taking a nap"
        import time
        time.sleep(60)
        
        print "stopping instance"
        gravity.stop()
        gravity2.stop()
        instance.stop()

    

    #run a single function job
    def test2(self):
        instance = self.start_nodes()
        def do_something(x, y):
            return x * y
        arguments = [8, 2]

        #job_id = instance.submit_function_job(self.do_something, arguments, node_label='VU')
        
        job_id = instance.submit_pickled_function_job("function pickle", "arguments pickle", node_label='VU')
        
        #result = instance.get_result(job_id)
        result = instance.get_pickled_function_job_result(job_id)
        
        instance.stop()
        
    #run a single script job
    def test3(self):
        instance = self.start_nodes()
        
        script_arguments = "--mass 3"

        job_id = instance.submit_script_job(script="some_script.py",
                                            arguments = script_arguments,
                                            script_dir = ".",
                                            node_label='VU')
        
        instance.stop()
        
    #run a parameter sweep
    def test4(self):
        instance = self.start_nodes()
        
        for model in ["small", "large", "huge"]:
            for mass in [1, 2, 3, 4, 5, 6]:
                input_file = model + ".particles"
                
                instance.submit_script_job(script="some_script.py",
                                           script_dir = ".",
                                           input_files=input_file,
                                           re_use_code_files = True,
                                           arguments = "--model " + input_file + " --mass " + str(mass))
                
        instance.wait_for_jobs()
        
        instance.stop()
        
