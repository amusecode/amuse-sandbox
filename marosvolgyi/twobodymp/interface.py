from amuse.community import *

class twobodympInterface(LegacyInterface):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        LegacyInterface.__init__(self, name_of_the_worker="twobodymp_worker", **keyword_arguments)
    
    @legacy_function
    def initialization():
        function = LegacyFunctionSpecification()  
        function.addParameter('precision', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def plot():
        function = LegacyFunctionSpecification()  
        function.addParameter('R', dtype='int32', direction=function.IN)
        function.addParameter('G', dtype='int32', direction=function.IN)
        function.addParameter('B', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def viewer():
        function = LegacyFunctionSpecification()  
        function.addParameter('view_', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_precision():
        function = LegacyFunctionSpecification()
        function.addParameter('precision', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_position():
        function = LegacyFunctionSpecification()
        function.addParameter('paritcle_id', dtype ='int32', direction=function.IN)
        function.addParameter('x_', dtype='float64', direction=function.IN)
        function.addParameter('y_', dtype='float64', direction=function.IN)
        function.addParameter('z_', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_position():
        function = LegacyFunctionSpecification()
        function.addParameter('paritcle_id', dtype ='int32', direction=function.IN)
        function.addParameter('x_', dtype='float64', direction=function.OUT)
        function.addParameter('y_', dtype='float64', direction=function.OUT)
        function.addParameter('z_', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_reduced_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('mu_', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('paritcle_id', dtype ='int32', direction=function.IN)
        function.addParameter('mass_', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def commit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_center_of_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('cmx', dtype='float64', direction=function.OUT)
        function.addParameter('cmy', dtype='float64', direction=function.OUT)
        function.addParameter('cmz', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_kinetic_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('id_', dtype='int32', direction=function.IN)
        function.addParameter('Ek', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_velocity():
        function = LegacyFunctionSpecification()
        function.addParameter('paritcle_id', dtype ='int32', direction=function.IN)
        function.addParameter('vx_', dtype='float64', direction=function.IN)
        function.addParameter('vy_', dtype='float64', direction=function.IN)
        function.addParameter('vz_', dtype='float64', direction=function.IN)
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

    @legacy_function
    def get_position_mp():
        function = LegacyFunctionSpecification()
        function.addParameter('paritcle_id', dtype ='int32', direction=function.IN)
        function.addParameter('xmp', dtype='string', direction=function.INOUT)
        function.addParameter('ymp', dtype='string', direction=function.INOUT)
        function.addParameter('zmp', dtype='string', direction=function.INOUT)
        function.addParameter('precision', dtype ='int32', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def get_velocity_mp():
        function = LegacyFunctionSpecification()
        function.addParameter('paritcle_id', dtype ='int32', direction=function.IN)
        function.addParameter('vxmp', dtype='string', direction=function.INOUT)
        function.addParameter('vymp', dtype='string', direction=function.INOUT)
        function.addParameter('vzmp', dtype='string', direction=function.INOUT)
        function.addParameter('precision', dtype ='int32', direction=function.IN)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function
    def set_velocity_mp():
        function = LegacyFunctionSpecification()
        function.addParameter('paritcle_id', dtype ='int32', direction=function.IN)
        function.addParameter('vxmp', dtype='string', direction=function.IN)
        function.addParameter('vymp', dtype='string', direction=function.IN)
        function.addParameter('vzmp', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def evolve_system():
        function = LegacyFunctionSpecification()
        function.addParameter('new_time', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_particle():
        """
        Define a new particle in the stellar dynamics code. The particle is initialized with the provided
        mass, radius, position and velocity. This function returns an index that can be used to refer
        to this particle.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT, description =
            """
            An index assigned to the newly created particle.
            This index is supposed to be a local index for the code
            (and not valid in other instances of the code or in other codes)
            """
            )

        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The mass of the particle")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The radius of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.result_type = 'int32'
        function.result_doc = """ 0 - OK
            particle was created and added to the model
        -1 - ERROR
            only one particle could can be created"""
        return function

    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

class twobodymp(InCodeComponentImplementation):

    def __init__(self):
        InCodeComponentImplementation.__init__(self,  twobodympInterface())
    
