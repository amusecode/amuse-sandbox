from amuse.community import *

class supportInterface(LegacyInterface):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        LegacyInterface.__init__(self, name_of_the_worker="support_worker", **keyword_arguments)
    
    @legacy_function
    def add():
        function = LegacyFunctionSpecification()  
        function.addParameter('term1', dtype='float64', direction=function.IN)
        function.addParameter('term2', dtype='float64', direction=function.IN)
        function.addParameter('sum', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

class support(CodeInterface):

    def __init__(self):
        CodeInterface.__init__(self,  supportInterface())
    
