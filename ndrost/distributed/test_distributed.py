from amuse.community import *
from amuse.test.amusetest import TestWithMPI
from amuse.lab import *

from interface import DistributedAmuseInterface
from interface import DistributedAmuse

class DistributedTests(TestWithMPI):
    
    def start_nodes(self):
        instance = DistributedAmuse(redirection='none')
        instance.initialize_code()

        instance.new_resource(name='DAS4-VU',
                              hostname="fs0.das4.cs.vu.nl",
                              username="niels",
                              scheduler_type="sge", 
                              amuse_dir="/home/niels/amuse-svn",
                              )
        
        instance.new_resource(name='LGM',
                              hostname="fs.lgm.liacs.nl", 
                              amuse_dir='/home/niels/amuse-svn',
                              scheduler_type="local", 
                              username="niels")
        
        instance.new_resource(name='LGM-4', 
                              hostname="node004", 
                              gateways="fs.lgm.liacs.nl",
                              scheduler_type="local", 
                              amuse_dir='/var/local/amuse',
                              username="niels")
    
        instance.new_reservation(resource_name='DAS4-VU', node_count=5, time= 2|units.hour, node_label='VU')
        
        instance.new_reservation(resource_name='LGM-4', node_count=1, time=2|units.hour, node_label='LGM')
    
        instance.wait_for_reservations()

        return instance

    def test1(self):
       instance = self.start_nodes()
    
       #gadget = Gadget(nr_of_workers=4, node_label='VU')
        
       #gadget2 = Gadget(nr_of_workers=4, nr_of_nodes=2, node_label='VU')

       # some interesting simulation using these workers

       instance.stop()

    def do_something(x, y):
        return x * y

    #run a single function job
    def test2(self):
        instance = self.start_nodes()
        
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
        
