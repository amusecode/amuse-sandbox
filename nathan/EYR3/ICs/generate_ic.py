import numpy
from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particles
from amuse.plot import plot, scatter, xlabel, ylabel, native_plot
#~from amuse.support.exceptions import AmuseException
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from amuse.ext.sph_to_star import convert_SPH_to_stellar_model
from amuse.community.kepler.interface import Kepler
from amuse.community.hop.interface import Hop

from amuse.ic.fractalcluster import new_fractal_cluster_model



    def relax(self, sph_model):
        gas_particles = sph_model.gas_particles
        hydro = self.new_hydrodynamics(gas_particles)
        hydro.parameters.artificial_viscosity_alpha = 0.0 # Viscous damping doesn't seem to be very important, but turned off just in case...
        channel_from_hydro = hydro.gas_particles.new_channel_to(gas_particles)
        channel_to_hydro = gas_particles.new_channel_to(hydro.gas_particles)
        few_dynamical_timescales = 3 * (2 * gas_particles.total_radius()**3 / (constants.G * gas_particles.total_mass())).sqrt()
        if self.verbose:
            print "Relaxing SPH model with {0} for {1}.".format(
                self.hydrodynamics.__name__, 
                few_dynamical_timescales.as_string_in(units.day))
        
        n_steps = 100
        if self.debug:
            monitor = dict(time=[]|units.day, kinetic=[]|units.J, potential=[]|units.J, thermal=[]|units.J)
        
#~        velocity_damp_factor = (1.0 - 2 * (few_dynamical_timescales / n_steps) / self.dynamical_timescale)
        velocity_damp_factor = 0.8
        for i_step, time in enumerate(few_dynamical_timescales * numpy.linspace(1.0/n_steps, 1.0, n_steps)):
            hydro.evolve_model(time)
            channel_from_hydro.copy_attributes(["mass","x","y","z","vx","vy","vz","u"])
            gas_particles.position -= gas_particles.center_of_mass()
#~            gas_particles.velocity = (gas_particles.velocity - gas_particles.center_of_mass_velocity()) * (i_step * 1.0 / n_steps)
            gas_particles.velocity = velocity_damp_factor * (gas_particles.velocity - gas_particles.center_of_mass_velocity())
            channel_to_hydro.copy_attributes(["mass","x","y","z","vx","vy","vz","u"])
            if self.debug:
                monitor["time"].append(time)
                monitor["kinetic"].append(hydro.kinetic_energy)
                monitor["potential"].append(hydro.potential_energy)
                monitor["thermal"].append(hydro.thermal_energy)
                
        hydro.stop()
        if self.debug:
            energy_evolution_plot(monitor["time"], monitor["kinetic"], monitor["potential"], monitor["thermal"])
        return gas_particles


if __name__ == "__main__":
    stars = new_fractal_cluster_model(N=1000, )

