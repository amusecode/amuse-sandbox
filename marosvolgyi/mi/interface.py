from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics

class MiInterface(CodeInterface,
                      LiteratureReferencesMixIn,
                      GravitationalDynamicsInterface,
                      StoppingConditionInterface):
    """
    N-body module with mixed 4th and 6th order Hermite integration scheme.


    """
    include_headers = ['worker_code.h', 'stopcond.h']

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="mi_worker",
                                 **options)
        LiteratureReferencesMixIn.__init__(self)

    def reinitialize_particles(self):
        self.recommit_particles()

    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_particle', dtype='int32',
                              direction=function.IN,
            description = "Index of particle to delete")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
           The particle was deleted
        """
        return function

#need to verify
    @legacy_function
    def set_smoothing():
        """
        Set the smoothing length basied on interaction type.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_particle1', dtype='int32',
                              direction=function.IN,
            description = "the index of the first particle")
        function.addParameter('index_of_particle2', dtype='int32',
                              direction=function.IN,
            description = "the index of the second particle")
        function.addParameter('smoothing_length', dtype='float64',
                              direction=function.IN,
            description = "set the smothing length between particle1 and particle2")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function

    @legacy_function
    def get_smoothing():
        """
        Get the smoothing length basied on interaction type.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_particle1', dtype='int32',
                              direction=function.IN,
            description = "the index of the first particle")
        function.addParameter('index_of_particle2', dtype='int32',
                              direction=function.IN,
            description = "the index of the second particle")
        function.addParameter('smoothing_length', dtype='float64',
                              direction=function.OUT,
            description = "set the smothing length between particle1 and particle2")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function

    @legacy_function
    def get_dt_dia():
        """
        Get the time interval between diagnostics output.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('dt_dia', dtype='float64',
                              direction=function.OUT,
            description = "time interval between diagnostic outputs")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_dt_dia():
        """
        Set the time interval between diagnostics output.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('dt_dia', dtype='float64',
                              direction=function.IN,
            description = "the time interval between diagnostics output")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_dt_param():
        """
        Get the timestep scaling factor.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('dt_dia', dtype='float64',
                              direction=function.OUT,
            description = "the timestep scaling factor")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_dt_param():
        """
        Set the timestep scaling factor.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('dt_dia', dtype='float64',
                              direction=function.IN,
            description = "the timestep scaling factor")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_time():
        """
        Get the current simulation time.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64',
                              direction=function.OUT,
            description = "the current simulation time")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_time():
        """
        Set the current simulation time.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64',
                              direction=function.IN,
            description = "the current simulation time")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_pair_detect_factor():
        """
        Get pair detection sphere radius factor (units particle radius).
        """
        function = LegacyFunctionSpecification()
        function.addParameter('pair_detect_factor', dtype='float64',
                              direction=function.OUT,
            description = "pair detection radius ")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def set_pair_detect_factor():
        """
        Set pair detection sphere radius factor (units particle radius).
        """
        function = LegacyFunctionSpecification()
        function.addParameter('pair_detect_factor', dtype='float64',
                              direction=function.IN,
            description = "pair detection radius ")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function




class Masaki(GravitationalDynamics):

##    __doc__ = MasakiDoc()

    def __init__(self, convert_nbody = None, **options):
        self.stopping_conditions = StoppingConditions(self)

        legacy_interface = MasakiInterface(**options)
        self.legacy_doc = legacy_interface.__doc__

        GravitationalDynamics.__init__(
            self,
            legacy_interface,
            convert_nbody,
            **options
        )

    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2",
            "set_eps2", 
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.0 | nbody_system.length * nbody_system.length
        )
#need to verify
#        object.add_method_parameter(
#            "get_smoothing",
#            "set_smoothing",
#            nbody_system.length * nbody_system.length, 
#            0.0 | nbody_system.length * nbody_system.length
#        )
        object.add_method_parameter(
            "get_dt_param",
            "set_dt_param",
            "dt_param",
            "timestep scaling factor", 
            units.none, 
            0.03 | units.none
        )
        object.add_method_parameter(
            "get_dt_dia",
            "set_dt_dia",
            "dt_dia", 
            "time interval between diagnostics output", 
            nbody_system.time,
            1.0 | nbody_system.time
        )
        object.add_method_parameter(
            "get_time",
            "set_time",
            "time",
            "current simulation time", 
            nbody_system.time, 
            0.0 | nbody_system.time
        )
        object.add_method_parameter(
            "get_pair_detect_factor",
            "set_pair_detect_factor",
            "pair_factor",
            "radius factor for pair detection", 
            units.none, 
            1.0 | units.none
        )

        self.stopping_conditions.define_parameters(object)

    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        object.add_method(
            'get_potential_at_point',
            (nbody_system.length, nbody_system.length,
             nbody_system.length, nbody_system.length),
            (nbody_system.potential, object.ERROR_CODE)
        )

#        object.add_method(
#            '',
#            (nbody_system.length, nbody_system.length,
#             nbody_system.length, nbody_system.length),
#            (nbody_system.potential, object.ERROR_CODE)
#        )

        object.add_method(
            'get_gravity_at_point',
            (nbody_system.length, nbody_system.length,
             nbody_system.length, nbody_system.length),
            (nbody_system.acceleration, nbody_system.acceleration,
             nbody_system.acceleration, object.ERROR_CODE)
        )
        
        self.stopping_conditions.define_methods(object)
        
    
    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)
        self.stopping_conditions.define_particle_set(object, 'particles')
