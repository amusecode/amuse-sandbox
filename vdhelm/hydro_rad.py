"""
   Minimalistic routine for running a radiative transfer code.
"""
from numpy import random

from amuse.units.optparse import OptionParser
from amuse.units import units, nbody_system
from amuse.datamodel.particles import Particle
from matplotlib import pyplot
from amuse import plot as aplot

from amuse.community.sphray.interface import SPHRay
from amuse.community.simplex.interface import SimpleXSplitSet
from amuse.community.fi.interface import Fi
from amuse.ext.molecular_cloud import ism_cube

from amuse.support.console import set_printing_strategy
set_printing_strategy("custom", preferred_units = [units.MSun, units.parsec, units.Myr])

def plot_3d(ism):
    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(ism.x.number, ism.y.number, ism.z.number, c=ism.xion, s=100)
    pyplot.show()
    pyplot.close()

def plot_subset(ism):
    subset = ism.select(lambda z : abs(z) < 1 | units.parsec, ["z"])
    aplot.scatter(subset.x, subset.y, c=subset.xion, s=100)
    pyplot.show()
    pyplot.close()

def output_results(ism):
    plot_subset(ism)

def create_source(Lstar):
    source=Particle()
    source.position = (0, 0, 0) | units.parsec
    source.mass = 1 | units.MSun
    source.flux = Lstar/(20. | units.eV)
    source.luminosity = source.flux
    source.rho = 1 | units.amu/units.cm**3
    source.xion = 0.0
    source.u = 81 | units.kms**2
    return source

def create_ism(N, boxsize):
    internal_energy = 81 | units.kms**2
    rho = 1 | units.amu/units.cm**3
    ism = ism_cube(N, boxsize/2., rho, internal_energy).result
    ism.flux = 0. | units.s**-1
    ism.xion = 0.0

    return ism

def create_codes(boxsize, end_time, total_mass):
    unit_converter = nbody_system.nbody_to_si(total_mass, end_time)
    hydro = Fi(unit_converter)
    radiative = SimpleXSplitSet()
    radiative.parameters.box_size=1.001*boxsize
    return hydro, radiative

def hydro_and_radiative_transfer(N, end_time, Lstar, boxsize):
    sync_timestep=0.001 | units.Myr
    plot_timestep=0.01 | units.Myr

    source = create_source(Lstar)
    ism = create_ism(N, boxsize)

    hydro, radiative = create_codes(boxsize, end_time, ism.mass.sum())

    hydro.gas_particles.add_particles(ism)

    channels = []
    channels.append(hydro.gas_particles.new_channel_to(ism))
    channels[0].copy()

    radiative.src_particles.add_particle(source)
    radiative.gas_particles.add_particles(ism)

    channels.append(radiative.gas_particles.new_channel_to(ism))
    channels.append(hydro.gas_particles.new_channel_to(radiative.gas_particles))
    radiative_feedback = radiative.gas_particles.new_channel_to(hydro.gas_particles, attributes=["u"])

    time = 0 | units.Myr
    while time <= end_time:
        for channel in channels:
            channel.copy()

        print "evolving to time", time
        hydro.evolve_model(time)
        radiative_feedback.copy()
        radiative.evolve_model(time)

        if time % plot_timestep < sync_timestep:
            output_results(ism)

        time += sync_timestep

def parse_arguments():
    parser = OptionParser()
    parser.add_option("-N", dest="N", type="int", default=200,
        help="The number of particles representing the gas [%default].")
    parser.add_option("-t", dest="end_time", unit=units.Myr, default=0.01, type='float',
        help="The end time of the simulation [%default %unit].")
    parser.add_option("-L", dest="Lstar", unit=units.LSun, default = 100, type='float',
        help="Luminosity of ionizing source [%default %unit].")
    parser.add_option("-s", dest="boxsize", unit=units.parsec, default=10, type='float',
        help="The total size of the simulation box [%default %unit].")

    options, args = parser.parse_args()
    return options.__dict__

if __name__ == "__main__":
    options = parse_arguments()
    hydro_and_radiative_transfer(**options)
