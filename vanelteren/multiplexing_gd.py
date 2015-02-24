from amuse.community.interface import gd


from amuse.units import nbody_system
from amuse.units import generic_unit_converter
from amuse.community.interface import common
from amuse.rfi.core import *

from amuse import datamodel

from amuse.community.hermite0.interface import Hermite

from amuse.community.mercury.interface import Mercury

import numpy

from collections import namedtuple

from amuse.units import units

class MultiplexingGravitationalDynamicsInterface(PythonCodeInterface, gd.GravitationalDynamicsInterface, gd.GravityFieldInterface):
    
    def __init__(self, **options):
        PythonCodeInterface.__init__(self, implementation_factory = MultiplexingGravitationalDynamicsImplementation, **options)
    
    @legacy_function
    def new_multiplexed_particle():
        """
        Define a new particle in the stellar dynamics code. The particle is initialized with the provided
        mass, radius, position and velocity. This function returns an index that can be used to refer
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

        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The radius of the particle", default = 0)
        function.addParameter('index_of_the_set', dtype='int32', direction=function.IN, description = "The index of the multiplexed set", default = 0)
        function.result_type = 'int32'
        function.result_doc = """ 0 - OK
            particle was created and added to the model
        -1 - ERROR
            particle could not be created"""
        return function



    @legacy_function
    def get_mass():
        """
            Retrieve the mass of a particle. Mass is a scalar property of a particle,
            this function has one OUT argument.
            """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
                              description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.OUT, description = "The current mass of the particle")
        function.result_type = 'int32'
        function.must_handle_array = True
        function.result_doc = """
            0 - OK
                particle was removed from the model
            -1 - ERROR
                particle could not be found
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
        function.result_type = 'int32'
        function.must_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was retreived
        -1 - ERROR
            particle could not be found
        """
        return function



    @legacy_function
    def get_position():
        """
        Retrieve the position vector of a particle. Position is a vector property,
        this function has 3 OUT arguments.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.result_type = 'int32'
        function.must_handle_array = True
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            not yet implemented
        """
        return function


    @legacy_function
    def get_velocity():
        """
        Retrieve the velocity vector of a particle. Position is a vector property,
        this function has 3 OUT arguments.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the velocity from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('vx', dtype='float64', direction=function.OUT, description = "The current x component of the position vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.OUT, description = "The current y component of the position vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.OUT, description = "The current z component of the position vector of the particle")
        function.result_type = 'int32'
        function.must_handle_array = True
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            not yet implemented
        """
        return function



    @legacy_function
    def get_state():
        """
        Retrieve the current state of a particle. The *minimal* information of a stellar
        dynamics particle (mass, radius, position and velocity) is returned.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.OUT, description = "The current mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('radius', dtype='float64', direction=function.OUT, description = "The current radius of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be found
        """
        return function


    @legacy_function
    def set_mass():
        """
        Update the mass of a particle. Mass is a scalar property of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
        function.result_type = 'int32'
        function.must_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        """
        return function


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
        function.must_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was retreived
        -1 - ERROR
            particle could not be found
        """
        return function


    @legacy_function
    def set_position():
        """
        Update the position of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.result_type = 'int32'
        function.must_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        """
        return function


    @legacy_function
    def set_velocity():
        """
        Set the velocity vector of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The current x component of the velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The current y component of the velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The current z component of the velocity vector of the particle")
        function.result_type = 'int32'
        function.must_handle_array = True
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            not yet implemented
        """
        return function



    @legacy_function
    def set_state():
        """
        Update the current state of a particle. The *minimal* information of a stellar
        dynamics particle (mass, radius, position and velocity) is updated.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The new radius of the particle", default = 0)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        -3 - ERROR
            not yet implemented
        """
        return function


    @legacy_function
    def delete_particle():
        """
        Remove the definition of particle from the code. After calling this function the particle is
        no longer part of the model evolution. It is up to the code if the index will be reused.
        This function is optional.
        """
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
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
        

    @legacy_function
    def set_code():
        """
        Update the name of the integration (multiplexed) code
        """
        function = LegacyFunctionSpecification()
        function.addParameter('code', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            code name was set
        -1 - ERROR
            code name is not known
        """
        return function


class SinglePointGravityFieldInterface(object):
    """
    Codes implementing the gravity field interface provide functions to
    calculate the force and potential energy fields at any point.
    """
    
    @legacy_function    
    def get_gravity_at_point():
        """
        Get the gravitational acceleration at the given points. To calculate the force on
        bodies at those points, multiply with the mass of the bodies
        """
        function = LegacyFunctionSpecification()  
        for x in ['eps','x','y','z']:
            function.addParameter(
              x, 
              dtype='float64', 
              direction=function.IN,
              unit=nbody_system.length
            )
        for x in ['ax','ay','az']:
            function.addParameter(
                x, 
                dtype='float64', 
                direction=function.OUT,     
                unit=nbody_system.acceleration
            )
        function.result_type = 'int32' 
        function.can_handle_array = True
        return function
        
    @legacy_function    
    def get_potential_at_point():
        """
        Determine the gravitational potential on any given point
        """
        function = LegacyFunctionSpecification()  
        for x in ['eps','x','y','z']:
            function.addParameter(
                x, 
                dtype='float64', 
                direction=function.IN,
                unit=nbody_system.length
            )
        for x in ['phi']:
            function.addParameter(
                x, 
                dtype='float64', 
                direction=function.OUT,
                unit=nbody_system.potential
            )
        function.result_type = 'int32'
        function.can_handle_array = True
        return function


class MultiplexingGravitationalDynamicsImplementation(object):
    
    def __init__(self):
        self.keys_generator = datamodel.BasicUniqueKeyGenerator(1)
        self.particles = datamodel.Particles()
        self.name_of_the_code = "Hermite"
        CodeDefinition = namedtuple('CodeDefinition', ['code', 'mass_unit', 'speed_unit', 'length_unit', 'time_unit'])
        self.code_names_to_code_classes = {
            "Hermite": CodeDefinition(Hermite,nbody_system.mass, nbody_system.speed, nbody_system.length, nbody_system.time),
            "Mercury": CodeDefinition(Mercury,units.MSun, units.AUd, units.AU, units.day)
        }
        self.code = None
        self.begin_time = 0
        self.time = 0
        self.mass_unit = nbody_system.mass
        self.speed_unit = nbody_system.speed
        self.length_unit = nbody_system.length
        self.time_unit = nbody_system.time
        

    def new_multiplexed_particle(self, index_of_the_particle, mass, x, y, z, vx, vy, vz, radius, index_of_the_set):
        new_particles = datamodel.Particles(len(mass), keys_generator = self.keys_generator)
        new_particles.mass = mass | self.mass_unit
        new_particles.x = x | self.length_unit
        new_particles.y = y | self.length_unit
        new_particles.z = z | self.length_unit
        new_particles.vx = vx | self.speed_unit
        new_particles.vy = vy | self.speed_unit
        new_particles.vz = vz | self.speed_unit
        new_particles.radius = radius | self.length_unit
        new_particles.index_of_the_set = index_of_the_set
        self.particles.add_particles(new_particles)
        index_of_the_particle.value = new_particles.key
        return 0
        


    def get_mass(self, index_of_the_particle, mass):
        subset = self.particles._subset(index_of_the_particle)
        if not len(subset) == len(index_of_the_particle):
            mass.value = [0] * len(index_of_the_particle)
            return -1
        else:
            mass.value = subset.mass.value_in(self.mass_unit)
            return 0
        


    def get_position(self, index_of_the_particle, x, y, z):
        subset = self.select_particles(index_of_the_particle)
        x.value = subset.x.value_in(self.length_unit)
        y.value = subset.y.value_in(self.length_unit)
        z.value = subset.z.value_in(self.length_unit)
        return 0
        


    def get_velocity(self, index_of_the_particle, vx, vy, vz):
        subset = self.select_particles(index_of_the_particle)
        vx.value = subset.vx.value_in(self.speed_unit)
        vy.value = subset.vy.value_in(self.speed_unit)
        vz.value = subset.vz.value_in(self.speed_unit)
        return 0
        


    def get_state(self, index_of_the_particle, mass, x, y, z, vx, vy, vz, radius):
        print index_of_the_particle
        subset = self.select_particles(index_of_the_particle)
        if not len(subset) == len(index_of_the_particle):
            return -1
        else:
            mass.value = subset.mass.value_in(self.mass_unit)
            x.value = subset.x.value_in(self.length_unit)
            y.value = subset.y.value_in(self.length_unit)
            z.value = subset.z.value_in(self.length_unit)
            vx.value = subset.vz.value_in(self.speed_unit)
            vy.value = subset.vy.value_in(self.speed_unit)
            vz.value = subset.vz.value_in(self.speed_unit)
            radius.value = subset.radius.value_in(self.length_unit)
            return 0
        

    def select_particles(self,keys):
        return self.particles._subset(keys)
        



    def set_mass(self, index_of_the_particle, mass):
        self.particles._subset(index_of_the_particle).mass = mass | self.mass_unit
        return 0
        



    def set_position(self, index_of_the_particle, x, y, z):
        subset = self.select_particles(index_of_the_particle)
        subset.x = x | self.length_unit
        subset.y = y | self.length_unit
        subset.z = z | self.length_unit
        return 0
        



    def set_velocity(self, index_of_the_particle, vx, vy, vz):
        subset = self.select_particles(index_of_the_particle)
        subset.vx = vx | self.speed_unit
        subset.vy = vy | self.speed_unit
        subset.vz = vz | self.speed_unit
        return 0
        



    def set_state(self, index_of_the_particle, mass, x, y, z, vx, vy, vz, radius):
        subset = self.select_particles(index_of_the_particle)
        subset.mass = mass | self.mass_unit
        subset.x = x | self.length_unit
        subset.y = y | self.length_unit
        subset.z = z | self.length_unit
        subset.vz = vz | self.speed_unit
        subset.vy = vy | self.speed_unit
        subset.vz = vz | self.speed_unit
        subset.radius = radius | self.length_unit
        return 0
        


    def delete_particle(self, index_of_the_particle):
        subset = self.select_particles(index_of_the_particle)
        self.particles.remove_particles(subset)
        return 0
        



    def initialize_code(self):
        self.keys_generator = datamodel.BasicUniqueKeyGenerator(1)
        self.particles = datamodel.Particles()
        return 0
        



    def cleanup_code(self):
        self.keys_generator = datamodel.BasicUniqueKeyGenerator(1)
        self.particles = datamodel.Particles()
        if not self.code is None:
            self.code.stop()
            self.code = None
        return 0
        



    def commit_particles(self):
        return 0
        



    def recommit_particles(self):
        return 0
        



    def commit_parameters(self):
        code_definition = self.code_names_to_code_classes[self.name_of_the_code]
        self.factory = code_definition.code
        
        self.mass_unit = code_definition.mass_unit
        self.speed_unit = code_definition.speed_unit
        self.length_unit = code_definition.length_unit
        self.time_unit = code_definition.time_unit
        return 0
        



    def recommit_parameters(self):
        return 0
        



    def set_begin_time(self, time):
        self.begin_time = time
        return 0
        



    def get_begin_time(self, time):
        time.value = self.begin_time
        



    def get_eps2(self, epsilon_squared):
        epsilon_squared.value = self.epsilon_squared
        return 0
        



    def set_eps2(self, epsilon_squared):
        self.epsilon_squared = epsilon_squared
        return 0
        



    def get_code(self, code):
        code.value = self.name_of_the_code
        return 0
        



    def set_code(self, code):
        self.name_of_the_code = code
        return 0
        



    def evolve_model(self, time):
        print "self.time_unit", self.time_unit
        if self.code is None:
            factory = self.code_names_to_code_classes[self.name_of_the_code].code
            self.code = factory(redirection="null")#debugger="ddd")
            self.code.parameters.epsilon_squared = self.epsilon_squared | self.length_unit**2
            
        set_indices = numpy.unique(self.particles.index_of_the_set)
        for index in set_indices:
            particles = self.particles[self.particles.index_of_the_set == index]
            print "index:", index, len(particles)
            self.code.reset()
            self.code.parameters.begin_time = self.time | self.time_unit
            self.code.particles.add_particles(particles)
            
            self.code.evolve_model(time | self.time_unit)
            print self.code.model_time
            particles.velocity = self.code.particles.velocity
            particles.position = self.code.particles.position
        self.time = time
        return 0
        



    def get_time(self, time):
        time.value = self.begin_time
        



    def synchronize_model(self):
        return 0
        



class MultiplexingGravitationalDynamicsCode(gd.GravitationalDynamics):
    

    def __init__(self, unit_converter = None,  **options):
        
        remote_code = MultiplexingGravitationalDynamicsInterface()
        gd.GravitationalDynamics.__init__(self, remote_code, unit_converter, **options)

    def define_properties(self, object):
        object.add_property('get_time', public_name = "model_time")
        

    def define_methods(self, object):
        gd.GravitationalDynamics.define_methods(self, object)
        
        object.add_method(
            "new_multiplexed_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
                object.INDEX,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )
        
        object.add_method(
            "get_eps2",
            (),
            (nbody_system.length * nbody_system.length, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_eps2",
            (nbody_system.length * nbody_system.length, ),
            (object.ERROR_CODE,)
        )


    def define_particle_sets(self, object):
        object.define_set('particles', 'index_of_the_particle')
        object.set_new('particles', 'new_multiplexed_particle')
        object.set_delete('particles', 'delete_particle')
        object.add_setter('particles', 'set_state')
        object.add_getter('particles', 'get_state')
        object.add_setter('particles', 'set_mass')
        object.add_getter('particles', 'get_mass', names = ('mass',))
        object.add_setter('particles', 'set_position')
        object.add_getter('particles', 'get_position')
        object.add_setter('particles', 'set_velocity')
        object.add_getter('particles', 'get_velocity')
        object.add_setter('particles', 'set_radius')
        object.add_getter('particles', 'get_radius')
        
        

    def define_state(self, object): 
        gd.GravitationalDynamics.define_state(self, object)               
        object.add_method('EDIT', 'new_multiplexed_particle')
        object.add_method('UPDATE', 'new_multiplexed_particle')
        object.add_transition('RUN', 'UPDATE', 'new_multiplexed_particle', False)
        


    def define_parameters(self, object):
        object.add_method_parameter(
            "get_eps2",
            "set_eps2", 
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            default_value = 0.0 | nbody_system.length * nbody_system.length
        )
        object.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )



class MultiplexingMercuryCode(MultiplexingGravitationalDynamicsCode):
    

    def __init__(self, unit_converter = None,  **options):
        MultiplexingGravitationalDynamicsCode.__init__(self, unit_converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(1|units.day, 1|units.AU, 1|units.MSun), **options)

