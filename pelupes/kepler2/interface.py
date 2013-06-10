from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics, GravityFieldCode
from amuse.support.options import option
from amuse.units import units
from amuse.datamodel import Particles
import os.path

class Kepler2Interface(CodeInterface,
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
                               name_of_the_worker = "kepler2_worker",
                               **options)

    @legacy_function
    def get_central_mass():
        """
        Return the current time of the system.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('central_mass', dtype='float64', direction=function.OUT,
                              unit = nbody_system.mass)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_central_mass():
        """
        Return the current time of the system.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('central_mass', dtype='float64', direction=function.IN,
                              unit = nbody_system.mass)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_semi_major_axis():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index', dtype='i', direction=function.IN)
        function.addParameter('semi_major_axis', dtype='float64', direction=function.OUT,
                              unit = nbody_system.length)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_eccentricity():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index', dtype='i', direction=function.IN)
        function.addParameter('eccentricity', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_specific_orbital_energy():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index', dtype='i', direction=function.IN)
        function.addParameter('specific_orbital_energy', dtype='float64', direction=function.OUT,
                              unit = nbody_system.specific_energy)
        function.result_type = 'int32'
        return function


    @legacy_function
    def get_next_radial_crossing_time():
        function = LegacyFunctionSpecification()
        function.can_handle_array=True
        function.addParameter('index', dtype='i', direction=function.IN)
        function.addParameter('target_radius', dtype='float64', direction=function.IN, unit=nbody_system.length)
        function.addParameter('next_radial_crossing_time', dtype='float64', direction=function.OUT, unit=nbody_system.time)
        function.result_type = 'int32'
        return function

    def get_gravity_at_point(self,radius,x,y,z):
        mass,err=self.get_central_mass()
        dr2=(x**2+y**2+z**2+radius**2)
        dr=dr2**0.5
        ax=-mass*x/(dr2*dr)
        ay=-mass*y/(dr2*dr)
        az=-mass*z/(dr2*dr)
        return ax,ay,az

    def get_potential_at_point(self,radius,x,y,z):
        mass,err=self.get_central_mass()
        dr2=(x**2+y**2+z**2+radius**2)
        dr=dr2**0.5
        phi=-mass/dr
        return phi

class Kepler2OrbitersOnly(GravitationalDynamics, GravityFieldCode):

    def __init__(self, unit_converter = None,  **options):
        self.unit_converter = unit_converter
        
        CommonCode.__init__(self,
                               Kepler2Interface(**options),
                               **options)

    def define_converter(self, object):
        if not self.unit_converter is None:
            object.set_converter(self.unit_converter.as_converter_from_si_to_generic())
    
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_central_mass", 
            "set_central_mass",
            "central_mass", 
            "central mass parameter", 
            default_value = 0 | nbody_system.mass
        )
         
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
        GravitationalDynamics.define_particle_sets(self, object)
        object.add_getter('particles', 'get_semi_major_axis')
        object.add_getter('particles', 'get_eccentricity')
        object.add_getter('particles', 'get_specific_orbital_energy')
        object.add_method('particles', 'get_next_radial_crossing_time')
  

class Kepler2(Kepler2OrbitersOnly):
    def __init__(self, *args, **kargs):
        Kepler2OrbitersOnly.__init__(self, *args, **kargs)
        self._particles=Particles()
        self.particles_accessed=True
        self.central_particle=None

    @property
    def orbiters_astro_centric(self):
        if self.particles_accessed:
            self.commit_particles()
        return self.overridden().particles

    @property
    def particles(self):
        if not self.particles_accessed:
            channel=self.orbiters_astro_centric.new_channel_to(self._particles)
            channel.copy_attributes(["x","y","z","vx","vy","vz","mass"])        
            orbiters=self.orbiters_astro_centric.get_intersecting_subset_in(self._particles)
            orbiters.position+=self.central_particle.position
            orbiters.velocity+=self.central_particle.velocity        
            self.particles_accessed=True
        return self._particles
        
    def commit_particles(self):
        self.particles_accessed=False

        N=len(self._particles)
        if N<1:
          raise Exception("too few particles")

        ic=numpy.argmax(self._particles.mass)        
        self.central_particle=self._particles[ic]
        
        self.parameters.central_mass=self.central_particle.mass
        
        orbiters=self._particles.copy()
        orbiters.remove_particle(self.central_particle)
                
        orbiters.position=orbiters.position-self.central_particle.position
        orbiters.velocity=orbiters.velocity-self.central_particle.velocity

        if (orbiters.mass.sum()/self.central_particle.mass) > 1.e-7:
          raise Exception("too heavy orbiters")
        
        orbiters.synchronize_to(self.orbiters_astro_centric)        
        channel=orbiters.new_channel_to(self.orbiters_astro_centric)
        channel.copy_attributes(["x","y","z","vx","vy","vz","mass"])

        return 0

    def recommit_particles(self):
        self.commit_particles()
                        
    def evolve_model(self, tend):
        if self.particles_accessed:
            self.recommit_particles()

        dt=tend-self.model_time
        
        self.overridden().evolve_model(tend)
        
        self.central_particle.position+=self.central_particle.velocity*dt
        self.particles_accessed=False

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
        
        
from amuse.couple.parallel_stellar_evolution import ParallelStellarEvolution
        
class ParallelKepler2(object):
    def __init__(self, **kargs):
        n=kargs.pop("number_of_workers")
        self.worker=ParallelStellarEvolution( Kepler2OrbitersOnly, n, **kargs)        
        self._particles=Particles()
        self.particles_accessed=True
        self.central_particle=None

    @property
    def parameters(self):
        return self.worker.parameters

    @property
    def model_time(self):
        return self.worker.model_time

    @property
    def orbiters_astro_centric(self):
        if self.particles_accessed:
            self.commit_particles()
        return self.worker.particles

    @property
    def particles(self):
        if not self.particles_accessed:
            channel=self.orbiters_astro_centric.new_channel_to(self._particles)
            channel.copy_attributes(["x","y","z","vx","vy","vz","mass"])        
            orbiters=self.orbiters_astro_centric.get_intersecting_subset_in(self._particles)
            orbiters.position+=self.central_particle.position
            orbiters.velocity+=self.central_particle.velocity        
            self.particles_accessed=True
        return self._particles
        
    def commit_particles(self):
        self.particles_accessed=False

        N=len(self._particles)
        if N<1:
          raise Exception("too few particles")

        ic=numpy.argmax(self._particles.mass)        
        self.central_particle=self._particles[ic]
        
        self.parameters.central_mass=self.central_particle.mass
        
        orbiters=self._particles.copy()
        orbiters.remove_particle(self.central_particle)
                
        orbiters.position=orbiters.position-self.central_particle.position
        orbiters.velocity=orbiters.velocity-self.central_particle.velocity

        if (orbiters.mass.sum()/self.central_particle.mass) > 1.e-7:
          raise Exception("too heavy orbiters")
        
        orbiters.synchronize_to(self.orbiters_astro_centric)        
        channel=orbiters.new_channel_to(self.orbiters_astro_centric)
        channel.copy_attributes(["x","y","z","vx","vy","vz","mass"])

    def recommit_particles(self):
        self.commit_particles()
                        
    def evolve_model(self, tend):
        if self.particles_accessed:
            self.recommit_particles()

        dt=tend-self.model_time
        
        self.worker.evolve_model(tend)
        
        self.central_particle.position+=self.central_particle.velocity*dt
        self.particles_accessed=False

    def get_gravity_at_point(self,radius,x,y,z):
        if self.particles_accessed:
            self.recommit_particles()
        xx=x-self.central_particle.x
        yy=y-self.central_particle.y
        zz=z-self.central_particle.z
        return self.worker.get_gravity_at_point(radius,xx,yy,zz)

    def get_potential_at_point(self,radius,x,y,z):
        if self.particles_accessed:
            self.recommit_particles()
        xx=x-self.central_particle.x
        yy=y-self.central_particle.y
        zz=z-self.central_particle.z
        return self.worker.get_potential_at_point(radius,xx,yy,zz)

