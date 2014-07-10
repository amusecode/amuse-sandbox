import numpy

from amuse.units import units, nbody_system
from amuse.datamodel.particles import Particles
from amuse.support.console import set_printing_strategy

from amuse.ic.plummer import new_plummer_model
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.ic.salpeter import new_salpeter_mass_distribution

from amuse.community.fi.interface import Fi
from amuse.community.huayno.interface import Huayno
from amuse.couple.bridge import Bridge

from matplotlib import pyplot
from amuse import plot as aplot

converter = nbody_system.nbody_to_si(100.|units.MSun, 2.|units.parsec)
cmapdict = {'red': ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)),
            'green': ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
            'blue': ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))}

pyplot.register_cmap("RedBlack", data=cmapdict)

star_positions = [
    # E
    [-4.5, 4.5],
    [-4.5, 4.3],
    [-4.5, 4.1],
    [-4.5, 3.9],
    [-4.5, 3.7],
    [-4.5, 3.5],
    [-4.5, 3.3],
    [-4.5, 3.1],

    [-4.3, 4.5],
    [-4.1, 4.5],
    [-3.9, 4.5],

    [-4.3, 3.1],
    [-4.1, 3.1],
    [-3.9, 3.1],

    [-4.3, 3.8],
    [-4.1, 3.8],

    # M
    [3.1, -3.7],
    [3.1, -3.9],
    [3.1, -4.1],
    [3.1, -4.3],
    [3.1, -4.5],

    [4.3, -3.7],
    [4.3, -3.9],
    [4.3, -4.1],
    [4.3, -4.3],
    [4.3, -4.5],

    [3.3, -3.9],
    [3.5, -4.1],
    [3.7, -4.3],
    [3.9, -4.1],
    [4.1, -3.9],

    # Heart
    [0.0, -3.0],
    [0.15, -2.8],
    [0.3, -2.6],
    [0.45, -2.4],
    [0.6, -2.2],
    [0.75, -2.0],
    [0.9, -1.8],
    [1.05, -1.6],
    [1.2, -1.4],
    [1.35, -1.2],
    [1.5, -1.0],
    [1.65, -0.8],
    [1.8, -0.6],
    [1.95, -0.4],
    [2.05, -0.2],
    [2.1, -0.0],
    [2.1, 0.2],
    [2.05, 0.4],
    [1.95, 0.6],
    [1.8, 0.8],
    [1.65, 0.95],
    [1.45, 1.05],
    [1.25, 1.1],
    [1.05, 1.1],
    [0.9, 1.05],
    [0.75, 0.95],
    [0.6, 0.8],
    [0.45, 0.65],
    [0.3, 0.5],
    [0.15, 0.35],
    [0.0, 0.2],

    [-0.15, -2.8],
    [-0.3, -2.6],
    [-0.45, -2.4],
    [-0.6, -2.2],
    [-0.75, -2.0],
    [-0.9, -1.8],
    [-1.05, -1.6],
    [-1.2, -1.4],
    [-1.35, -1.2],
    [-1.5, -1.0],
    [-1.65, -0.8],
    [-1.8, -0.6],
    [-1.95, -0.4],
    [-2.05, -0.2],
    [-2.1, -0.0],
    [-2.1, 0.2],
    [-2.05, 0.4],
    [-1.95, 0.6],
    [-1.8, 0.8],
    [-1.65, 0.95],
    [-1.45, 1.05],
    [-1.25, 1.1],
    [-1.05, 1.1],
    [-0.9, 1.05],
    [-0.75, 0.95],
    [-0.6, 0.8],
    [-0.45, 0.65],
    [-0.3, 0.5],
    [-0.15, 0.35],

    # Arrow
    [-3.2, 2.2],
    [-3.0, 2.0],
    [-2.8, 1.8],
    [-2.6, 1.6],
    [-2.4, 1.4],
    [-2.2, 1.2],
    [-2.0, 1.0],
    [-1.8, 0.8],
    [-1.6, 0.6],
    [-1.4, 0.4],
    [-1.2, 0.2],
    [-1.0, 0.0],
    [-0.8, -0.2],
    [-0.6, -0.4],
    [-0.4, -0.6],

    [1.0, -2.0],
    [1.2, -2.2],
    [1.4, -2.4],
    [1.6, -2.6],
    [1.8, -2.8],
    [2.0, -3.0],
    [2.2, -3.2],
    [2.4, -3.4],
    [2.6, -3.6],

    # Feather
    [-3.2, 2.4],
    [-3.2, 2.6],
    [-3.4, 2.2],
    [-3.6, 2.2],

    [-2.8, 2.0],
    [-2.8, 2.2],
    [-3.0, 1.8],
    [-3.2, 1.8],

    [-3.0, 2.4],
    [-3.4, 2.0],

    # Point
    [2.6, -3.4],
    [2.6, -3.2],
    [2.6, -3.0],
    [2.4, -3.6],
    [2.2, -3.6],
    [2.0, -3.6],

    ] | units.parsec

def create_stars():
    number = len(star_positions)
    stars = new_plummer_model(number, converter)
    stars.mass = new_salpeter_mass_distribution(number) #TODO set numpy random seed
    stars.x = star_positions[:,0]
    stars.y = star_positions[:,1]
    # stars.move_to_center()
    return stars

def create_gas():
    gas = new_plummer_gas_model(500, converter)
    # gas.position =
    gas.move_to_center()
    return gas

def setup_codes(stars, gas):
    grav = Huayno(converter)
    stars = grav.particles.add_particles(stars)
    hydro = Fi(converter)
    gas = hydro.gas_particles.add_particles(gas)

    bridge = Bridge(use_threading=False)
    bridge.add_system(hydro, (grav, ) )
    bridge.add_system(grav)

    hydro.parameters.timestep = 0.1|units.Myr
    bridge.timestep = 0.2|units.Myr
    grav.parameters.epsilon_squared = (0.01|units.parsec)**2

    return bridge, stars, gas

def plot(stars, gas, i):
    aplot.pynbody_column_density_plot(gas, width=10|units.parsec, cmap="RedBlack", vmin=20, vmax=22)
    aplot.scatter(stars.x, stars.y, color='yellow', s=20)

    aplot.xlim([-5, 5]|units.parsec)
    aplot.ylim([-5, 5]|units.parsec)
    aplot.xlabel("x")
    aplot.ylabel("y")
    pyplot.savefig("m_plots/plot_{:=04}0.png".format(i))

def evolve(end_time=100|units.Myr, dt = 0.2|units.Myr):
    bridge, stars, gas = setup_codes(create_stars(), create_gas())

    t = 0.|units.Myr
    i = int(end_time / dt)
    while t < end_time:
        print "evolving to {} (i={})".format(t, i)
        bridge.evolve_model(t)

        plot(stars, gas, i)

        t += dt
        i -= 1

if __name__ == '__main__':
    set_printing_strategy("custom", preferred_units = [units.MSun, units.parsec, units.Myr])
    evolve()

