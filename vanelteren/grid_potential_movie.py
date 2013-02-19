from amuse.lab import *

from matplotlib import pyplot

def hydro_grid_in_potential_well(mass = 1 | units.MSun, length = 100 | units.AU):
    converter = nbody_system.nbody_to_si(mass, length)
    
    # calculate density in field based on solar wind
    # gives a very low number
    molar_mass_hydrogen_proton = 1 | units.g / units.mol
    density_hydrogen_in_stellar_wind = 10 | 1 / units.cm**3
    particles_per_mol =  6.022e23 | 1 / units.mol
    density_hydrogen_in_stellar_wind_in_moles = density_hydrogen_in_stellar_wind / particles_per_mol
    density_gas = 100 * (density_hydrogen_in_stellar_wind_in_moles * molar_mass_hydrogen_proton).as_quantity_in(units.MSun / units.AU**3)
    
    
    # override with higher number for plotting
    #density_gas = 1e-3 | units.MSun / units.AU**3
    
    instance=Athena(converter)
    instance.initialize_code()
    instance.parameters.nx = 500
    instance.parameters.ny = 500
    instance.parameters.nz = 1
    instance.parameters.length_x = length 
    instance.parameters.length_y = length
    instance.parameters.length_z = length
    instance.parameters.x_boundary_conditions = ("outflow","outflow")
    instance.parameters.y_boundary_conditions = ("outflow","outflow")
    instance.parameters.z_boundary_conditions = ("outflow","outflow")
    
    
    instance.set_has_external_gravitational_potential(1)
    
    instance.commit_parameters()
    
    grid_in_memory = instance.grid.copy()
    grid_in_memory.rho  = density_gas
    pressure =  1e-12 | units.Pa
    R = 8.314  | units.J * units.K**-1 * units.mol **-1
    temperature = (pressure / density_gas) * molar_mass_hydrogen_proton / R
    
    print "temperature of the interplanetary gas:", temperature.as_quantity_in(units.K)
    
    grid_in_memory.energy =  pressure / (instance.parameters.gamma- 1)
    channel = grid_in_memory.new_channel_to(instance.grid)
    channel.copy()
    
    instance.initialize_grid()
    particle = Particle(
        mass = mass,
        position = length * [0.5, 0.5, 0.5],
        velocity = [0.0, 0.0, 0.0] | units.kms
    )
    
    gravity = Hermite(converter)
    dx = (grid_in_memory.x[1][0][0] - grid_in_memory.x[0][0][0]).as_quantity_in(units.AU)
    gravity.parameters.epsilon_squared = (10 * dx)**2 
    gravity.particles.add_particle(particle)
    
    
    potential = gravity.get_potential_at_point(
        0 * instance.potential_grid.x.flatten(),
        instance.potential_grid.x.flatten(),
        instance.potential_grid.y.flatten(),
        instance.potential_grid.z.flatten()
    )
    
    potential =  potential.reshape(instance.potential_grid.x.shape)
    instance.potential_grid.potential =  potential
    
    dt = 1 | units.yr
    t = 0 | units.yr
    while t < 200 | units.yr:
        t += dt
        print t.value_in(units.yr), instance.get_timestep().value_in(units.yr)
        instance.evolve_model(t)
        print t.value_in(units.yr), instance.get_timestep().value_in(units.yr)
        value_to_plot = instance.grid.rho[...,...,0].value_in(units.MSun / units.AU**3)
        #value_to_plot = potential[...,...,0].value_in(potential.unit)
        plot_grid(value_to_plot, t.value_in(units.yr))
    
def plot_grid(x, t):
    figure = pyplot.figure(figsize=(6,6))
    plot = figure.add_subplot(1,1,1)
    mappable = plot.imshow(x, origin = 'lower')
    pyplot.colorbar(mappable)
    figure.savefig('grid_potential_{0:03}.png'.format(int(t)))
    
    
if __name__ == '__main__':
    hydro_grid_in_potential_well()
        
