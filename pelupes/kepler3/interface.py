from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics, GravityFieldCode
from amuse.support.options import option
from amuse.units import units
from amuse.datamodel import Particles
import os.path

class KeplerInterface(CodeInterface,
                       GravitationalDynamicsInterface):
    """
    Kepler orbit manipulation functions, imported from Starlab.
    Initialize an orbit from mass, pos, and vel, or mass, semi-major
    axis and eccentricity, and allow the user to manipulate the
    resulting structure.  Most Starlab functionality is currently
    exposed.
    """

    # Interface specification.

    include_headers = ['interface.h']
    
    def __init__(self, **options):
        CodeInterface.__init__(self,
                               name_of_the_worker = "kepler_worker",
                               **options)

    @legacy_function
    def get_central_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='i', direction=function.IN, default=0)        
        function.addParameter('mass', dtype='float64', direction=function.OUT,
                              unit = nbody_system.mass)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_central_mass():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('mass', dtype='float64', direction=function.IN,
                              unit = nbody_system.mass)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_central_pos():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='i', direction=function.IN, default=0)        
        function.addParameter('x', dtype='float64', direction=function.OUT,
                              unit = nbody_system.length)
        function.addParameter('y', dtype='float64', direction=function.OUT,
                              unit = nbody_system.length)
        function.addParameter('z', dtype='float64', direction=function.OUT,
                              unit = nbody_system.length)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_central_pos():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.addParameter('y', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.addParameter('z', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_central_vel():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='i', direction=function.IN, default=0)        
        function.addParameter('vx', dtype='float64', direction=function.OUT,
                              unit = nbody_system.speed)
        function.addParameter('vy', dtype='float64', direction=function.OUT,
                              unit = nbody_system.speed)
        function.addParameter('vz', dtype='float64', direction=function.OUT,
                              unit = nbody_system.speed)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_central_vel():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='i', direction=function.IN)
        function.addParameter('vx', dtype='float64', direction=function.IN,
                              unit = nbody_system.speed)
        function.addParameter('vy', dtype='float64', direction=function.IN,
                              unit = nbody_system.speed)
        function.addParameter('vz', dtype='float64', direction=function.IN,
                              unit = nbody_system.speed)
        function.result_type = 'int32'
        return function

    @legacy_function
    def new_central_particle():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='i', direction=function.OUT)        
        function.addParameter('mass', dtype='float64', direction=function.IN,
                              unit = nbody_system.mass)
        function.addParameter('x', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.addParameter('y', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.addParameter('z', dtype='float64', direction=function.IN,
                              unit = nbody_system.length)
        function.addParameter('vx', dtype='float64', direction=function.IN,
                              unit = nbody_system.speed)
        function.addParameter('vy', dtype='float64', direction=function.IN,
                              unit = nbody_system.speed)
        function.addParameter('vz', dtype='float64', direction=function.IN,
                              unit = nbody_system.speed)
        function.addParameter('radius', dtype='float64', direction=function.IN,
                              unit = nbody_system.length,default=0.)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_central_particle():
        function = LegacyFunctionSpecification()
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'int32'
        return function



    @legacy_function
    def get_semi_major_axis():
        function = LegacyFunctionSpecification()
        function.must_handle_array=True
        function.addParameter('index', dtype='i', direction=function.IN)
        function.addParameter('semi_major_axis', dtype='float64', direction=function.OUT,
                              unit = nbody_system.length)
        function.addParameter('number_of_points', 'i', function.LENGTH)           
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eccentricity():
        function = LegacyFunctionSpecification()
        function.must_handle_array=True
        function.addParameter('index', dtype='i', direction=function.IN)
        function.addParameter('eccentricity', dtype='float64', direction=function.OUT)
        function.addParameter('number_of_points', 'i', function.LENGTH)           
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_specific_orbital_energy():
        function = LegacyFunctionSpecification()
        function.must_handle_array=True
        function.addParameter('index', dtype='i', direction=function.IN)
        function.addParameter('specific_orbital_energy', dtype='float64', direction=function.OUT,
                              unit = nbody_system.specific_energy)
        function.addParameter('number_of_points', 'i', function.LENGTH)           
        function.result_type = 'int32'
        return function


    @legacy_function
    def get_next_radial_crossing_time():
        function = LegacyFunctionSpecification()
        function.must_handle_array=True
        function.addParameter('index', dtype='i', direction=function.IN)
        function.addParameter('target_radius', dtype='float64', direction=function.IN, unit=nbody_system.length)
        function.addParameter('next_radial_crossing_time', dtype='float64', direction=function.OUT, unit=nbody_system.time)
        function.addParameter('number_of_points', 'i', function.LENGTH)           
        function.result_type = 'int32'
        return function

    def get_gravity_at_point(self,radius,x,y,z):
        mass,err=self.get_central_mass()
        xc,yc,zc,err=self.get_central_pos()
        dr2=((x-xc)**2+(y-yc)**2+(z-zc)**2+radius**2)
        dr=dr2**0.5
        ax=-mass*(x-xc)/(dr2*dr)
        ay=-mass*(y-yc)/(dr2*dr)
        az=-mass*(z-zc)/(dr2*dr)
        return ax,ay,az

    def get_potential_at_point(self,radius,x,y,z):
        mass,err=self.get_central_mass()
        xc,yc,zc,err=self.get_central_pos()
        dr2=((x-xc)**2+(y-yc)**2+(z-zc)**2+radius**2)
        dr=dr2**0.5
        phi=-mass/dr
        return phi

class Kepler(GravitationalDynamics, GravityFieldCode):

    def __init__(self, unit_converter = None,  **options):
        self.unit_converter = unit_converter
        
        CommonCode.__init__(self,
                               KeplerInterface(**options),
                               **options)

    def define_converter(self, object):
        if not self.unit_converter is None:
            object.set_converter(self.unit_converter.as_converter_from_si_to_generic())
    
         
    def define_methods(self, object):
        GravitationalDynamics.define_methods(self, object)
        object.add_method(
            'get_gravity_at_point',
            (
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
            ),
            (
                nbody_system.acceleration,
                nbody_system.acceleration,
                nbody_system.acceleration,
            )
        )
        object.add_method(
            'get_potential_at_point',
            (
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
            ),
            (
                nbody_system.potential,
            )
        )

      
    def define_particle_sets(self, object):
        object.define_super_set('particles', ['central_particle','orbiters'], 
            index_to_default_set = 1)

        object.define_set('central_particle', 'index_of_the_particle')
        object.set_new('central_particle', 'new_central_particle')
        object.set_delete('central_particle', 'delete_central_particle')
        object.add_setter('central_particle', 'set_central_mass')
        object.add_getter('central_particle', 'get_central_mass')
        object.add_setter('central_particle', 'set_central_pos')
        object.add_getter('central_particle', 'get_central_pos')
        object.add_setter('central_particle', 'set_central_vel')
        object.add_getter('central_particle', 'get_central_vel')

        object.define_set('orbiters', 'index_of_the_particle')
        object.set_new('orbiters', 'new_particle')
        object.set_delete('orbiters', 'delete_particle')
        object.add_setter('orbiters', 'set_state')
        object.add_getter('orbiters', 'get_state')
        object.add_setter('orbiters', 'set_mass')
        object.add_getter('orbiters', 'get_mass')
        object.add_setter('orbiters', 'set_position')
        object.add_getter('orbiters', 'get_position')
        object.add_setter('orbiters', 'set_velocity')
        object.add_getter('orbiters', 'get_velocity')
        object.add_getter('orbiters', 'get_semi_major_axis')
        object.add_getter('orbiters', 'get_eccentricity')
        object.add_getter('orbiters', 'get_specific_orbital_energy')
        object.add_method('orbiters', 'get_next_radial_crossing_time')
  

    def get_gravity_at_point(self,radius,x,y,z):
        if self.particles_accessed:
            self.recommit_particles()
        xx=x-self.central_particle.x
        yy=y-self.central_particle.y
        zz=z-self.central_particle.z
        return self.overridden().get_gravity_at_point(radius,xx,yy,zz)

    def get_potential_at_point(self,radius,x,y,z):
        if self.particles_accessed:
            self.recommit_particles()
        xx=x-self.central_particle.x
        yy=y-self.central_particle.y
        zz=z-self.central_particle.z
        return self.overridden().get_potential_at_point(radius,xx,yy,zz)
