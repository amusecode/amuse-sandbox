"""
Plummer 2d model generator

This module contains a function used to create Plummer (1911) models, which 
follow a spherically symmetric density profile of the form:
rho = c * (1 + r**2)**(-5/2)
"""

import numpy
import numpy.random

from math import pi, sqrt

from amuse.units import nbody_system
from amuse import datamodel
from amuse.ic import plummer

__all__ = ["new_plummer_model_2d"]

class MakePlummerModel2D(plummer.MakePlummerModel):
   
    def new_positions_spherical_coordinates(self):
        pi2 = pi * 2
        radius = self.calculate_radius_uniform_distribution()
        phi = self.random.uniform(0.0,pi2, (self.number_of_particles,1))
        return (radius,phi)

    def new_velocities_spherical_coordinates(self, radius):
        pi2 = pi * 2
        x,y = self.new_xy_for_velocity()
        velocity = x * sqrt(2.0) * numpy.power( 1.0 + radius*radius, -0.25)
        phi = self.random.uniform(0.0,pi2, (self.number_of_particles,1))
        return (velocity,phi)

    def coordinates_from_spherical(self, radius, phi):
        x = radius * numpy.cos( phi )
        y = radius * numpy.sin( phi )
        z = radius * 0
        return (x,y,z)

    def new_model(self):
        m = numpy.zeros((self.number_of_particles,1)) + (1.0 / self.number_of_particles)
        radius, phi = self.new_positions_spherical_coordinates()
        position =  numpy.hstack(self.coordinates_from_spherical(radius, phi))
        radius, phi = self.new_velocities_spherical_coordinates(radius)
        velocity = numpy.hstack(self.coordinates_from_spherical(radius , phi))
        position = position / 1.695
        velocity = velocity / sqrt(1 / 1.695)
        return (m, position, velocity)

    @property
    def result(self):
        masses = numpy.ones(self.number_of_particles) / self.number_of_particles
        radius, phi = self.new_positions_spherical_coordinates()
        x,y,z =  self.coordinates_from_spherical(radius, phi)
        radius, phi = self.new_velocities_spherical_coordinates(radius)
        vx,vy,vz = self.coordinates_from_spherical(radius, phi)
 
        result = datamodel.Particles(self.number_of_particles)
        result.mass = nbody_system.mass.new_quantity(masses)
        result.x = nbody_system.length.new_quantity(x.reshape(self.number_of_particles)/1.695)
        result.y = nbody_system.length.new_quantity(y.reshape(self.number_of_particles)/1.695)
        result.z = nbody_system.length.new_quantity(z.reshape(self.number_of_particles)/1.695)
        result.vx = nbody_system.speed.new_quantity(vx.reshape(self.number_of_particles) / sqrt(1/1.695))
        result.vy = nbody_system.speed.new_quantity(vy.reshape(self.number_of_particles) / sqrt(1/1.695))
        result.vz = nbody_system.speed.new_quantity(vz.reshape(self.number_of_particles) / sqrt(1/1.695))
        result.radius = 0 | nbody_system.length

        result.move_to_center()
        if self.do_scale:
            result.scale_to_standard()

        if not self.convert_nbody is None:
            result = datamodel.ParticlesWithUnitsConverted(result, self.convert_nbody.as_converter_from_si_to_generic())
            result = result.copy()

        return result

def new_plummer_model_2D(number_of_particles, *list_arguments, **keyword_arguments):
    """
    Create a plummer sphere with the given number of particles. Returns
    a set of stars with equal mass and positions and velocities distributed
    to fit a plummer star distribution model. The model is centered around the
    origin. Positions and velocities are optionally scaled such that the kinetic and
    potential energies are 0.25 and -0.5 in nbody-units, respectively.

    :argument number_of_particles: Number of particles to include in the plummer sphere
    :argument convert_nbody:  When given will convert the resulting set to SI units
    :argument radius_cutoff: Cutoff value for the radius (defaults to 22.8042468)
    :argument mass_cutoff: Mass percentage inside radius of 1
    :argument do_scale: scale the result to exact nbody units (M=1, K=0.25, U=-0.5)
    """
    uc = MakePlummerModel2D(number_of_particles, *list_arguments, **keyword_arguments)
    return uc.result

