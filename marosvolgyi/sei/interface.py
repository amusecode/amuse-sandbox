from amuse.community import *

class SeiInterface(CodeInterface):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="sei_worker", **keyword_arguments)
    
    @legacy_function
    def initialization():
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.addParameter('vx', dtype='float64', direction=function.IN)
        function.addParameter('vy', dtype='float64', direction=function.IN)
        function.addParameter('vz', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.addParameter('x', dtype='float64', direction=function.OUT)
        function.addParameter('y', dtype='float64', direction=function.OUT)
        function.addParameter('z', dtype='float64', direction=function.OUT)
        function.addParameter('vx', dtype='float64', direction=function.OUT)
        function.addParameter('vy', dtype='float64', direction=function.OUT)
        function.addParameter('vz', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def evolve():
        function = LegacyFunctionSpecification()
        function.addParameter('t_end', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_velocity():
        function = LegacyFunctionSpecification()
        function.addParameter('paritcle_id', dtype ='int32', direction=function.IN)
        function.addParameter('vx_', dtype='float64', direction=function.OUT)
        function.addParameter('vy_', dtype='float64', direction=function.OUT)
        function.addParameter('vz_', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

class twobodymp(InCodeComponentImplementation):

    def __init__(self):
        InCodeComponentImplementation.__init__(self,  twobodympInterface())
    
