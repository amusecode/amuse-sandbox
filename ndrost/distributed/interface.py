import threading
import sys
import logging
import socket


from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface
from amuse.community.interface.common import CommonCode
from amuse.units import units
from amuse.support import options

logger = logging.getLogger(__name__)

class OutputHandler(threading.Thread):
    
    def __init__(self, stream, port):
        threading.Thread.__init__(self)
        self.stream = stream

        logging.getLogger("channel").debug("output handler connecting to distributed code")
        
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        
        address = ('localhost', port)
        
        try:
            self.socket.connect(address)
        except:
            raise exceptions.CodeException("Could not connect to distributed code at " + str(address))
        
        self.socket.setsockopt(socket.SOL_TCP, socket.TCP_NODELAY, 1)
        
        self.socket.sendall('TYPE_OUTPUT'.encode('utf-8'))

        #fetch ID of this connection
        
        result = SocketMessage()
        result.receive(self.socket)
        
        self.id = result.strings[0]
        
        self.daemon = True
        self.start()
        
    def run(self):
        
        while True:
            logging.getLogger("channel").debug("receiving data for output")
            data = self.socket.recv(1024)
            
            if len(data) == 0:
                logging.getLogger("channel").debug("end of output", len(data))
                return
            
            logging.getLogger("channel").debug("got %d bytes", len(data))
            
            self.stream.write(data)

class DistributedAmuseInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn):
    """
	Distributed Amuse Code
    
        .. [#] The Distributed Amuse project is a collaboration between Sterrewacht Leiden and The Netherlands eScience Center.
    """

    classpath = ['.', 'worker.jar', 'src/dist/*']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="distributed_worker_java", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
        
        port = self.get_worker_port()
        
        #logging.basicConfig(level=logging.DEBUG)
        
        logger.debug("running on port %d", port)

#        self.stdoutHandler = OutputHandler(sys.stdout, port)
#        self.stderrHandler = OutputHandler(sys.stderr, port)

        options.GlobalOptions.instance().override_value_for_option("channel_type", "ibis")
        options.GlobalOptions.instance().override_value_for_option("port", port)


    @option(choices=['mpi','remote','ibis', 'sockets'], sections=("channel",))
    def channel_type(self):
        return 'sockets'
    
    @legacy_function
    def get_worker_port():
        """
        Returns the server socket port of the code. Used by the distributed channel
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def new_resource():
        """
        Define a new resource. This function returns an index that can be used to refer
        to this resource.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_resource', dtype='int32', direction=function.OUT)
        function.addParameter("name", dtype='string', direction=function.IN)
        function.addParameter("location", dtype='string', direction=function.IN)
        function.addParameter("amuse_dir", dtype='string', direction=function.IN)
        function.addParameter("scheduler_type", dtype='string', direction=function.IN, default="fork")
        function.addParameter('start_hub', dtype='int32', direction=function.IN, default=-1)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_resource():
        """
        Remove the definition of resource from the code. After calling this function the resource is
        no longer part of the model evolution. It is up to the code if the index will be reused.
        This function is optional.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_resource', dtype='int32', direction=function.IN,
            description = "Index of the resource to be removed. This index must have been returned by an earlier call to :meth:`new_resource`")

        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            resource was removed from the model
        -1 - ERROR
            resource could not be removed
        -2 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def new_reservation():
        """
        Reserve one or more nodes for later use by the simulation.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('reservation_id', dtype='int32', direction=function.OUT)
        function.addParameter("resource_name", dtype='string', direction=function.IN)
        function.addParameter("queue_name", dtype='string', direction=function.IN, default="")
        function.addParameter("node_count", dtype='int32', direction=function.IN, default = 1)
        function.addParameter("time", dtype='int32', direction=function.IN, unit = units.minute, default = 60)
        function.addParameter("slots", dtype='int32', direction=function.IN, default = 1)
        function.addParameter("node_label", dtype='string', direction=function.IN, default = "default")
        function.addParameter('count', dtype='int32', direction=function.LENGTH)

        function.result_type = 'int32'
        return function

        
    @legacy_function
    def delete_reservation():
        """
        Delete a reservation.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('reservation_id', dtype='int32', direction=function.IN)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    
    @legacy_function
    def wait_for_reservations():
        """
        Wait until all reservations are started, and all nodes are available to run jobs and/or workers
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
        
    @legacy_function
    def submit_pickled_function_job():
        """
        Submit a job, specified by a pickle of the function, and a pickle of the arguments.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('job_id', dtype='int32', direction=function.OUT)
        function.addParameter('function', dtype='string', direction=function.IN)
        function.addParameter('arguments', dtype='string', direction=function.IN)
        function.addParameter("node_label", dtype='string', direction=function.IN, default = "default")
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_pickled_function_job_result():
        """
        Get a result of a picked function job. Will block until the result is available
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('job_id', dtype='int32', direction=function.IN)
        function.addParameter('result', dtype='string', direction=function.OUT)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def submit_script_job():
        """
        Submit a job, specified by a script
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('job_id', dtype='int32', direction=function.OUT)
        function.addParameter('script', dtype='string', direction=function.IN)
        function.addParameter('arguments', dtype='string', direction=function.IN)
        function.addParameter('script_dir', dtype='string', direction=function.IN)
        function.addParameter("node_label", dtype='string', direction=function.IN, default = "default")
        function.addParameter("re_use_code_files", dtype='int32', direction=function.IN, default = 0)
        function.addParameter('count', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def wait_for_jobs():
        """
        Wait until all jobs are done.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    def cleanup_code(self):
        del options.GlobalOptions.instance().overriden_options["channel_type"]
        return 0
    
    
    
class DistributedAmuse(CommonCode):

    def __init__(self, **options):
        CommonCode.__init__(self,  DistributedAmuseInterface(**options), **options)
    
        
#    def store_view(self, description=""):
#        self.overridden().store_view(str(description))
#        
#    def define_parameters(self, object):
#        object.add_boolean_parameter(
#            "get_use_star_shader_flag",
#            "set_use_star_shader_flag",
#            "use_star_shader",
#            "Use-star-shader flag. False means: plain spheres.",
#            True
#        )
#        object.add_boolean_parameter(
#            "get_use_octree_for_gas_flag",
#            "set_use_octree_for_gas_flag",
#            "use_octree_for_gas",
#            "Use-octree-for-gas flag. True means: gas resources are divided over "
#                "octree cells, and these cells will be visualized instead.",
#            False
#        )
    
