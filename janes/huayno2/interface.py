from amuse.community import *
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics
from amuse.support.units import units

# http://stackoverflow.com/questions/36932/whats-the-best-way-to-implement-an-enum-in-python
class inttypes:
    EVOLVE_SHARED=1
    EVOLVE_PASS_KDK=2
    EVOLVE_HOLD_KDK=3
    EVOLVE_BRIDGE_KDK=4
    EVOLVE_EXTRAPOLATE=5
    EVOLVE_PASS_DKD=7
    EVOLVE_HOLD_DKD=8
    EVOLVE_HOLD_DKD_CC=9
    EVOLVE_SHARED4=10
    EVOLVE_HOLD_DKD_CC4=11
    EVOLVE_OK2=12
    EVOLVE_OK4=13


class HuaynoInterface(CodeInterface,GravitationalDynamicsInterface):
    include_headers = ['worker_code.h']

    EVOLVE_SHARED=1
    EVOLVE_EXTRAPOLATE=5
    EVOLVE_PASS_KDK=2
    EVOLVE_PASS_DKD=7
    EVOLVE_HOLD_KDK=3
    EVOLVE_HOLD_DKD=8
    EVOLVE_BRIDGE_KDK=4
    EVOLVE_HOLD_DKD_CC=9
    EVOLVE_SHARED4=10
    EVOLVE_HOLD_DKD_CC4=11
    
    MODE_OPENCL='opencl'
    MODE_OPENMP='openmp'
        
    def name_of_worker(self,mode):
        if mode==self.MODE_OPENCL:
            return 'huayno_worker_cl'
        if mode==self.MODE_OPENMP:
            return 'huayno_worker_mp'
        return 'huayno_worker'
      
    def __init__(self, mode=None, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_worker(mode), **options) 
        
    @legacy_function    
    def new_particle():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.OUT)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_time():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function      
    def commit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_kinetic_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('kinetic_energy', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_potential_energy():
        function = LegacyFunctionSpecification()
        function.addParameter('potential_energy', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def initialize_code():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_state():
        function = LegacyFunctionSpecification()   
        function.can_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['mass','radius','x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def evolve():
        function = LegacyFunctionSpecification()  
        function.addParameter('time_end', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_timestep_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('time_param', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_timestep_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('time_param', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_inttype_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('inttype', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_inttype_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('inttype', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_eps2_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('eps2', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def set_eps2_parameter():
        function = LegacyFunctionSpecification()
        function.addParameter('eps2', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    def set_eps2(self, e):
        return self.set_eps2_parameter(e)

    def get_eps2(self):
        return self.get_eps2_parameter()
    
    @legacy_function    
    def get_evolve_statistics():
        function = LegacyFunctionSpecification()   
        function.addParameter('ttot', dtype='i', direction=function.OUT)
        function.addParameter('ktot', dtype='i', direction=function.OUT)
        function.addParameter('dtot', dtype='i', direction=function.OUT)
        function.addParameter('tstot', dtype='i', direction=function.OUT)
        function.addParameter('kstot', dtype='i', direction=function.OUT)
        function.addParameter('dstot', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function      
    def get_sys_initial_timestep():
        function = LegacyFunctionSpecification()
        function.addParameter('timestep', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def get_interaction_fw_timestep():
        function = LegacyFunctionSpecification()   
        function.addParameter('id_i', dtype='i', direction=function.IN)
        function.addParameter('id_j', dtype='i', direction=function.IN)
        function.addParameter('timestep', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function

    @legacy_function    
    def reverse_velocities():
        function = LegacyFunctionSpecification()   
        function.result_type = 'i'
        return function

    @legacy_function
    def get_conserved_quantities():
        function = LegacyFunctionSpecification()   
        function.addParameter('e', dtype='d', direction=function.OUT)
        function.addParameter('px', dtype='d', direction=function.OUT)
        function.addParameter('py', dtype='d', direction=function.OUT)
        function.addParameter('pz', dtype='d', direction=function.OUT)
        function.addParameter('Lx', dtype='d', direction=function.OUT)
        function.addParameter('Ly', dtype='d', direction=function.OUT)
        function.addParameter('Lz', dtype='d', direction=function.OUT)
        function.addParameter('cx', dtype='d', direction=function.OUT)
        function.addParameter('cy', dtype='d', direction=function.OUT)
        function.addParameter('cz', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
        
class Huayno(GravitationalDynamics):

    EVOLVE_SHARED=1
    EVOLVE_EXTRAPOLATE=5
    EVOLVE_PASS_KDK=2
    EVOLVE_PASS_DKD=7
    EVOLVE_HOLD_KDK=3
    EVOLVE_HOLD_DKD=8
    EVOLVE_BRIDGE_KDK=4
    EVOLVE_HOLD_DKD_CC=9
    EVOLVE_SHARED4=10
    EVOLVE_HOLD_DKD_CC4=11

    def __init__(self, convert_nbody = None, **options):
        legacy_interface = HuaynoInterface(**options)
#        self.legacy_doc = legacy_interface.__doc__

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

        object.add_method_parameter(
            "get_timestep_parameter",
            "set_timestep_parameter", 
            "timestep_parameter", 
            "timestep parameter for gravity calculations", 
            units.none, 
            0.03 | units.none
        )

        object.add_method_parameter(
            "get_inttype_parameter",
            "set_inttype_parameter", 
            "inttype_parameter", 
            "integrator method to use", 
            units.none, 
            8 | units.none
        )

        object.add_method_parameter(
            "get_time",
            "set_time",
            "time",
            "current simulation time", 
            nbody_system.time, 
            0.0 | nbody_system.time
        )


    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)

        object.add_method("get_evolve_statistics", (),
            (units.none, units.none, units.none, units.none, units.none, units.none, object.ERROR_CODE)
        )

    def define_particle_sets(self, object):
        GravitationalDynamics.define_particle_sets(self, object)

    def define_state(self, object):
        GravitationalDynamics.define_state(self, object)
        object.add_method('PARAMETER_CHANGE_B', 'set_eps2_parameter')
        object.add_method('RUN', 'get_kinetic_energy')
        object.add_method('RUN', 'get_potential_energy')
