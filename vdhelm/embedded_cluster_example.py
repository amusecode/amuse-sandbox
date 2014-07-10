"""
    Evolve a small embedded cluster of massive stars with stellar wind.
    TODO: currently far too slow for an example
"""
from amuse.units.optparse import OptionParser
from amuse.units import units, quantities, nbody_system
from amuse.support.console import set_printing_strategy

from amuse.ic.plummer import new_plummer_model
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ext.stellar_wind import StarsWithSimpleWind as StarsWithWind
# from amuse.ext.stellar_wind import Stars_With_Accelerating_Wind as Stars_With_Wind

from amuse.community.seba.interface import SeBa
from amuse.community.fi.interface import Fi
from amuse.community.ph4.interface import ph4
from amuse.couple.bridge import Bridge

def plot_or_save(stars, gas, time):
    print "    ", time, "->", len(gas)

def evolve_embedded_cluster(N_stars, N_init_gas, radius, gas_part_mass, t_end, start_age, n_steps):
    converter = nbody_system.nbody_to_si(15*N_stars|units.MSun, radius)
    stars = new_plummer_model(N_stars, converter)
    stars.mass = new_salpeter_mass_distribution(N_stars, mass_min=15|units.MSun)

    gas_converter = nbody_system.nbody_to_si(N_init_gas*gas_part_mass, radius)
    gas = new_plummer_gas_model(N_init_gas, gas_converter)

    stellar = SeBa()
    stellar.particles.add_particles(stars)
    gravity = ph4(converter)
    stars = gravity.particles.add_particles(stars)
    hydro = Fi(gas_converter)
    gas = hydro.gas_particles.add_particles(gas)

    stellar_to_grav = stellar.particles.new_channel_to(gravity.particles, attributes=["mass", "radius"])
    stellar.evolve_model(start_age)
    stellar_to_grav.copy()

    wind = StarsWithWind(stars, gas_part_mass, evolving_stars=stellar.particles)
    bridge = Bridge(use_threading=False)
    bridge.add_system(hydro, (gravity, wind) )
    bridge.add_system(gravity, (hydro,) )
    bridge.add_system(wind)

    hydro.parameters.timestep = t_end / n_steps / 8.
    bridge.timestep = 2*hydro.parameters.timestep
    gravity.parameters.epsilon_squared = (10|units.RSun)**2


    for time in quantities.linspace(0.|units.yr, t_end, n_steps + 1):
        print "Evolving to time", time
        bridge.evolve_model(time)
        stellar.evolve_model(start_age + time)

        stellar_to_grav.copy()

        if wind.has_new_wind_particles():
            wind_gas = wind.create_wind_particles()
            gas.add_particles(wind_gas)

        plot_or_save(stars, gas, time)

def parse_arguments():
    parser = OptionParser()
    parser.add_option("-N", dest="N_stars", type="int", default = 3,
                      help="number of stars [%default]")
    parser.add_option("-n", dest="N_init_gas", type="int", default = 10,
                      help="number of initial gas particles [%default]")
    parser.add_option("-R", dest="radius",
                      unit=units.AU, type="float", default = 1|units.parsec,
                      help="Radius of the cluster [%default]")
    parser.add_option("-M", dest="gas_part_mass",
                      unit=units.MSun, type="float", default = 1e-11|units.MSun,
                      help="Mass of the gas particles [%default]")
    parser.add_option("-t", dest="t_end",
                      unit=units.yr, type="float", default = 10|units.yr,
                      help="end time of the simulation [%default]")
    parser.add_option("-a", dest="start_age",
                      unit=units.yr, type="float", default = 10|units.Myr,
                      help="initial age of the stars [%default]")
    parser.add_option("-s", dest="n_steps", type="int", default = 1000,
                      help="number of diagnostics time steps [%default]")
    options, arguments = parser.parse_args()

    return options.__dict__

if __name__ == '__main__':
    set_printing_strategy("custom", preferred_units = [units.MSun, units.parsec, units.yr])
    options = parse_arguments()
    evolve_embedded_cluster(**options)
