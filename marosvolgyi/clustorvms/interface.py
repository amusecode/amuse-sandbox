from amuse.community import *

class ClustorvmsInterface(CodeInterface):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="clustorvms_worker", **keyword_arguments)
    
    @legacy_function
    def echo_int():
        function = LegacyFunctionSpecification()  
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.addParameter('int_out', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function
    
    @legacy_function
    def set_numBits():
        function = LegacyFunctionSpecification()  
        function.addParameter('int_in', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
        
class Clustorvms(InCodeComponentImplementation):

    def __init__(self):
        InCodeComponentImplementation.__init__(self,  ClustorvmsInterface())
    
