from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface
from amuse.community.interface.common import CommonCode
from amuse.units import units

class DistributedAmuseInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn):
    """
    AstroTray is a 3D visualization package for AMUSE simulations.
    
        .. [#] The AstroTray 3D visualization project is a collaboration between Sterrewacht Leiden and The Netherlands eScience Center.
    """

    imports = ['nl.esciencecenter.estars.Code']
    classpath = '.:log4j.properties:lib/*:external/*:../*:external/jogl/*'
    cwd = 'src'
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="estars_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)

    @option(choices=['mpi','remote','ibis', 'sockets'], sections=("channel",))
    def channel_type(self):
        return 'sockets'
    
    @legacy_function
    def new_resource():
        """
        Define a new resource in the visualisation code. The particle is initialized with the provided
        radius, position and color. This function returns an index that can be used to refer
        to this particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_resource', dtype='int32', direction=function.OUT)
        function.addParameter("name", dtype='string', direction=function.IN)
        function.addParameter("uri", dtype='string', direction=function.IN)
#~        for par in ["x", "y", "z"]:
#~            function.addParameter(par, dtype='float64', unit=generic_unit_system.length, direction=function.IN, 
#~                description = "The initial position vector of the particle")
#~        function.addParameter('radius', dtype='float64', unit=generic_unit_system.length, direction=function.IN, description = "The radius of the particle")
#~        for par in ["red", "green", "blue"]:
#~            function.addParameter(par, dtype='float64', direction=function.IN, 
#~                description = "The RGB color of the particle")
#~        function.addParameter("alpha", dtype='float64', direction=function.IN, description = "The opacity of the particle", default = 1.0)
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_resource():
        """
        Remove the definition of particle from the code. After calling this function the particle is
        no longer part of the model evolution. It is up to the code if the index will be reused.
        This function is optional.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        #function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to be removed. This index must have been returned by an earlier call to :meth:`new_particle`")

        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be removed
        -2 - ERROR
            not yet implemented
        """
        return function



    @legacy_function
    def get_radius():
        """
        Retrieve the radius of a particle. Radius is a scalar property of a particle,
        this function has one OUT argument.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the radius of. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('radius', dtype='float64', direction=function.OUT, description = "The current radius of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        #function.can_handle_array = True
        function.must_handle_array = True
        return function


    @legacy_function
    def set_radius():
        """
        Set the radius of a particle. Radius is a scalar property of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the radius of. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('radius', dtype='float64', unit=generic_unit_system.length, direction=function.IN, description = "The new radius of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        #function.can_handle_array = True
        function.must_handle_array = True
        return function

    @legacy_function
    def commit_particles():
        """
        Let the code perform initialization actions
        after all particles have been loaded in the model.
        Should be called before the first evolve call and
        after the last new_particle call.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Model is initialized and evolution can start
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention
        """
        return function

    @legacy_function
    def recommit_particles():
        """
        Let the code perform initialization actions
        after all particles have been loaded in the model.
        Should be called before the first evolve call and
        after the last new_particle call.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Model is initialized and evolution can start
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention
        """
        return function

    
    @legacy_function
    def get_position():
        """
        Retrieve the position vector of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for par in ["x", "y", "z"]:
            function.addParameter(par, dtype='float64', unit=generic_unit_system.length, direction=function.OUT, 
                description = "The current position vector of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def set_position():
        """
        Update the position of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for par in ["x", "y", "z"]:
            function.addParameter(par, dtype='float64', unit=generic_unit_system.length, direction=function.IN, 
                description = "The new position vector of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def get_color():
        """
        Retrieve the RGB color vector of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for par in ["red", "green", "blue"]:
            function.addParameter(par, dtype='float64', direction=function.OUT, 
                description = "The current RGB color vector of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def set_color():
        """
        Update the RGB color of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for par in ["red", "green", "blue"]:
            function.addParameter(par, dtype='float64', direction=function.IN, 
                description = "The new RGB color vector of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def get_opacity():
        """
        Retrieve the alpha (opacity) of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter("alpha", dtype='float64', direction=function.OUT, 
            description = "The opacity of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def set_opacity():
        """
        Update the alpha (opacity) of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter("alpha", dtype='float64', direction=function.IN, 
            description = "The new opacity of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def get_use_star_shader_flag():
        function = LegacyFunctionSpecification()
        function.addParameter("use_star_shader_flag", dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_use_star_shader_flag():
        function = LegacyFunctionSpecification()
        function.addParameter("use_star_shader_flag", dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_use_octree_for_gas_flag():
        function = LegacyFunctionSpecification()
        function.addParameter("use_octree_for_gas_flag", dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_use_octree_for_gas_flag():
        function = LegacyFunctionSpecification()
        function.addParameter("use_octree_for_gas_flag", dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def store_view():
        """
        Store and view the current model, corresponding to the given description.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='string', direction=function.IN,
            description = "The description of the scene.")
        function.result_type = 'int32'
        return function
    
class DistributedAmuse(CommonCode):

    def __init__(self, **options):
        CommonCode.__init__(self,  DistributedAmuseInterface(**options), **options)
        
    def store_view(self, description=""):
        self.overridden().store_view(str(description))

    def define_particle_sets(self, object):
        object.define_set('resources', 'index_of_the_particle')
        object.set_new('resources', 'new_resource')
        object.set_delete('resources', 'delete_resource')
        
        object.add_getter('resources', 'get_name')
        object.add_setter('resources', 'set_name')
        object.add_getter('resources', 'get_color')
        object.add_setter('resources', 'set_color')
        object.add_getter('resources', 'get_opacity')
        object.add_setter('resources', 'set_opacity')
        object.add_getter('resources', 'get_radius')
        object.add_setter('resources', 'set_radius')

    def define_state(self, object): 
        CommonCode.define_state(self, object)   
        object.add_transition('END', 'INITIALIZED', 'initialize_code', False)    
        
        object.add_transition('INITIALIZED','EDIT','commit_parameters')
        object.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        object.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        object.add_transition('UPDATE','CHANGE_PARAMETERS_UPDATE','before_set_parameter', False)
        object.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        object.add_transition('CHANGE_PARAMETERS_UPDATE','UPDATE','recommit_parameters')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        object.add_method('CHANGE_PARAMETERS_UPDATE','before_set_parameter')
        
        object.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        object.add_method('CHANGE_PARAMETERS_UPDATE','before_get_parameter')
        object.add_method('RUN', 'before_get_parameter')
        object.add_method('EDIT', 'before_get_parameter')
        object.add_method('UPDATE','before_get_parameter')
        
        object.add_method('EDIT', 'new_star_particle')
        object.add_method('EDIT', 'new_gas_particle')
        object.add_method('EDIT', 'new_sphere_particle')
        object.add_method('EDIT', 'new_marker_particle')
        object.add_method('EDIT', 'delete_particle')
        object.add_method('UPDATE', 'new_star_particle')
        object.add_method('UPDATE', 'new_gas_particle')
        object.add_method('UPDATE', 'new_sphere_particle')
        object.add_method('UPDATE', 'new_marker_particle')
        object.add_method('UPDATE', 'delete_particle')
        object.add_transition('EDIT', 'RUN', 'commit_particles')
        object.add_transition('RUN', 'UPDATE', 'new_star_particle', False)
        object.add_transition('RUN', 'UPDATE', 'new_gas_particle', False)
        object.add_transition('RUN', 'UPDATE', 'new_sphere_particle', False)
        object.add_transition('RUN', 'UPDATE', 'new_marker_particle', False)
        object.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        object.add_transition('UPDATE', 'RUN', 'recommit_particles')
        
        object.add_method('RUN', 'store_view')
        object.add_method('RUN', 'set_position')
        object.add_method('RUN', 'set_color')
        object.add_method('RUN', 'set_opacity')
        object.add_method('RUN', 'get_position')
        object.add_method('RUN', 'get_color')
        object.add_method('RUN', 'get_opacity')
        
    def define_parameters(self, object):
        object.add_boolean_parameter(
            "get_use_star_shader_flag",
            "set_use_star_shader_flag",
            "use_star_shader",
            "Use-star-shader flag. False means: plain spheres.",
            True
        )
        object.add_boolean_parameter(
            "get_use_octree_for_gas_flag",
            "set_use_octree_for_gas_flag",
            "use_octree_for_gas",
            "Use-octree-for-gas flag. True means: gas particles are divided over "
                "octree cells, and these cells will be visualized instead.",
            False
        )
    
