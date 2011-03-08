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
    
