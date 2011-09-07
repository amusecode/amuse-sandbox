from amuse.community import *
from amuse.support.literature import LiteratureReferencesMixIn
class PACOInterface(CodeInterface, LiteratureReferencesMixIn):
    """
    PACO - Pattern AutoCOrrelation orbit classification scheme
           construct repeated patterns in cartesian crossings in
           position and velocity space

       .. [#] Faber, N., Flitti, F. Boily, C.M., Collete, C., 
       .. [#] Patsis, P.A., Portegies Zwart, SF., 2010, 
       .. [#] *Preprint submitted to MNRAS*
    """
    include_headers = ['worker_code.h', 'src/PACO.h']
    
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="worker_code", **options)
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def PartialPatternConstruct():
        """
        Construct patterns based on Cartesian axis changes.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('x0', 'd', function.IN)
        function.addParameter('y0', 'd', function.IN)
        function.addParameter('vx0', 'd', function.IN)
        function.addParameter('vy0', 'd', function.IN)
        function.addParameter('x1', 'd', function.IN)
        function.addParameter('y1', 'd', function.IN)
        function.addParameter('sign0', 'i', function.OUT)
        function.addParameter('sign1', 'i', function.OUT)
        return function



    def PatternConstruct(self, x, y, vx, vy):
        pattern = []
        for i in range(len(x)-1) :
            sign0, sign1 = self.PartialPatternConstruct(
                x[i], y[i], vx[i], vy[i], x[i+1],  y[i+1]) 
            if sign0 >= 0 :
                pattern.append('Y' if sign0 else 'X')
        if sign1>=0:
            pattern.append('Y' if sign1 else 'X')
        return "".join(pattern)

    def auto_correlate(self, pattern, treshold) : 
        import numpy as np
        size = len(pattern)
        data = np.zeros(size)
        for ip in range(size) :
            for id in range(size) :
                if pattern[id] == pattern[(ip+id)%size] :
                    data[ip] += 1
            data[ip] = data[ip]/size
        for id in range(len(data)):
            if data[id] >= treshold:
                data[id] = 1
            else :
                data[id] = 0
        return data

    @legacy_function
    def cleanup_buffers():
        function = LegacyFunctionSpecification()
        function.can_handle_array = False
        return function
    
    @legacy_function
    def determine_pattern():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('x', 'd', function.IN)
        function.addParameter('y', 'd', function.IN)
        function.addParameter('vx', 'd', function.IN)
        function.addParameter('vy', 'd', function.IN)
        function.addParameter('number_of_points', 'i', function.LENGTH)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_upsilon():
        function = LegacyFunctionSpecification()
        function.addParameter('value', 'float32', function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_upsilon():
        function = LegacyFunctionSpecification()
        function.addParameter('value', 'float32', function.IN)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_psratio():
        function = LegacyFunctionSpecification()
        function.addParameter('value', 'float32', function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_repeating_pattern():
        function = LegacyFunctionSpecification()
        function.addParameter('value', 'string', function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_pattern_buffer():
        function = LegacyFunctionSpecification()
        function.addParameter('value', 'string', function.OUT)
        function.result_type = 'int32'
        return function 
        
    @legacy_function
    def get_delta_length():
        function = LegacyFunctionSpecification()
        function.addParameter('value', 'int32', function.OUT)
        function.result_type = 'int32'
        return function 
        
    @legacy_function
    def get_length_of_the_pattern():
        function = LegacyFunctionSpecification()
        function.addParameter('value', 'int32', function.OUT)
        function.result_type = 'int32'
        return function 

    @legacy_function
    def get_delta():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('output', 'int32', function.OUT)
        function.addParameter('value', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function 
        
    @legacy_function
    def get_autocorrelate_data():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('output', 'float32', function.OUT)
        function.addParameter('value', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function 
        
class PACO(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  PACOInterface(**options),  **options)
