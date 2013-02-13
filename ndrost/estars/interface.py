from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface
from amuse.community.interface.common import CommonCode
from amuse.units import units

length_unit = units.parsec
time_unit = units.Myr
type_unit = units.none

class eStarsInterface(CodeInterface, CommonCodeInterface):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="estars_worker", **keyword_arguments)
    
    @legacy_function
    def new_particle():
        """
        Define a new particle in the visualisation code. The particle is initialized with the provided
        type, position and color. This function returns an index that can be used to refer
        to this particle.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT, description =
            """
            An index assigned to the newly created particle.
            This index is supposed to be a local index for the code
            (and not valid in other instances of the code or in other codes)
            """
            )

        function.addParameter('type', dtype='int32', unit=type_unit, direction=function.IN, 
            description = "The type (gas, star, etc.) of the particle")
        for par in ["x", "y", "z"]:
            function.addParameter(par, dtype='float64', unit=length_unit, direction=function.IN, 
                description = "The initial position vector of the particle")
        function.addParameter('radius', dtype='float64', unit=length_unit, direction=function.IN, description = "The radius of the particle")
        for par in ["alpha", "red", "green", "blue"]:
            function.addParameter(par, dtype='float64', direction=function.IN, 
                description = "The ARGB color of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_particle():
        """
        Remove the definition of particle from the code. After calling this function the particle is
        no longer part of the model evolution. It is up to the code if the index will be reused.
        This function is optional.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to be removed. This index must have been returned by an earlier call to :meth:`new_particle`")
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



#    @legacy_function
#    def get_radius():
#        """
#        Retrieve the radius of a particle. Radius is a scalar property of a particle,
#        this function has one OUT argument.
#        """
#        function = LegacyFunctionSpecification()
#        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
#            description = "Index of the particle to get the radius of. This index must have been returned by an earlier call to :meth:`new_particle`")
#        function.addParameter('radius', dtype='float64', direction=function.OUT, description = "The current radius of the particle")
#        function.result_type = 'int32'
#        function.can_handle_array = True
#        function.result_doc = """
#        0 - OK
#            particle was found in the model and the information was retreived
#        -1 - ERROR
#            particle could not be found
#        """
#        return function


    @legacy_function
    def set_radius():
        """
        Set the radius of a particle. Radius is a scalar property of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the radius of. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The new radius of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was retreived
        -1 - ERROR
            particle could not be found
        """
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

    
#~    @legacy_function
#~    def get_position():
#~        """
#~        Retrieve the position vector of a particle.
#~        """
#~        function = LegacyFunctionSpecification()
#~        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
#~        for par in ["x", "y", "z"]:
#~            function.addParameter(par, dtype='float64', unit=length_unit, direction=function.OUT, 
#~                description = "The current position vector of the particle")
#~        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
#~        function.result_type = 'int32'
#~        function.must_handle_array = True
#~        return function
    
    @legacy_function
    def set_position():
        """
        Update the position of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for par in ["x", "y", "z"]:
            function.addParameter(par, dtype='float64', unit=length_unit, direction=function.IN, 
                description = "The new position vector of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
#~    @legacy_function
#~    def get_color():
#~        """
#~        Retrieve the RGB color vector of a particle.
#~        """
#~        function = LegacyFunctionSpecification()
#~        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
#~        for par in ["red", "green", "blue"]:
#~            function.addParameter(par, dtype='float64', direction=function.OUT, 
#~                description = "The current RGB color vector of the particle")
#~        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
#~        function.result_type = 'int32'
#~        function.must_handle_array = True
#~        return function
    
    @legacy_function
    def set_color():
        """
        Update the RGB color of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for par in ["red", "green", "blue", "alpha"]:
            function.addParameter(par, dtype='float64', direction=function.IN, 
                description = "The new RGBA color vector of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
#~    @legacy_function
#~    def get_type():
#~        """
#~        Retrieve the type of a particle.
#~        """
#~        function = LegacyFunctionSpecification()
#~        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
#~        function.addParameter('type', dtype='float64', unit=type_unit, direction=function.OUT, 
#~            description = "The type (gas, star, etc.) of the particle")
#~        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
#~        function.result_type = 'int32'
#~        function.must_handle_array = True
#~        return function
    
    @legacy_function
    def set_type():
        """
        Update the type of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('type', dtype='int32', unit=type_unit, direction=function.IN, 
            description = "The type (gas, star, etc.) of the particle")
        function.addParameter('npoints', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function
    
    @legacy_function
    def store_view():
        """
        Store and view the current model, corresponding to the given time.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', unit=time_unit, direction=function.IN,
            description = "The current model time.")
        function.result_type = 'int32'
        return function
    
class eStars(CommonCode):

    def __init__(self, **options):
        CommonCode.__init__(self,  eStarsInterface(**options), **options)
    

    def define_particle_sets(self, object):
        object.define_set('particles', 'index_of_the_particle')
        object.set_new('particles', 'new_particle')
#        object.set_delete('particles', 'delete_particle')
#~        object.add_getter('particles', 'get_type')
        object.add_setter('particles', 'set_type')
#~        object.add_getter('particles', 'get_position')
        object.add_setter('particles', 'set_position')
#~        object.add_getter('particles', 'get_color')
        object.add_setter('particles', 'set_color')
#        object.add_setter('particles', 'set_radius')
        object.add_setter('particles', 'set_radius')

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
#        object.add_method('EVOLVED','before_get_parameter')
        
        
        object.add_method('EDIT', 'new_particle')
#        object.add_method('EDIT', 'delete_particle')
        object.add_method('UPDATE', 'new_particle')
#        object.add_method('UPDATE', 'delete_particle')
        object.add_transition('EDIT', 'RUN', 'commit_particles')
        object.add_transition('RUN', 'UPDATE', 'new_particle', False)
 #       object.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        object.add_transition('UPDATE', 'RUN', 'recommit_particles')
        
#        object.add_transition('RUN', 'EVOLVED', 'store_view', False)
#        object.add_method('EVOLVED', 'store_view')
        object.add_method('RUN', 'store_view')
#        object.add_transition('EVOLVED','RUN', 'set_type')
#        object.add_transition('EVOLVED','RUN', 'set_position')
#        object.add_transition('EVOLVED','RUN', 'set_color')
        object.add_method('RUN', 'set_type')
        object.add_method('RUN', 'set_position')
        object.add_method('RUN', 'set_color')
#~        object.add_method('RUN', 'get_type')
#~        object.add_method('RUN', 'get_position')
#~        object.add_method('RUN', 'get_color')
        
