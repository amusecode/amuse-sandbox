import numpy
import pdb
from amuse.community.mercury.interface import MercuryWayWard
from amuse.ext.solarsystem import Solarsystem
from amuse.support.units import units
from amuse.support.units.values import VectorQuantity as Vq
from amuse.plot import *

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

def plotdata(planets):
    for planet in planets:
        t, x = planet.get_timeline_of_attribute_as_vector("x")
        t, y = planet.get_timeline_of_attribute_as_vector("y")
        plot(x, y, '.')
    native_plot.show()

def integrate_and_store():
    sun, planets = Solarsystem.new_solarsystem()
    #timerange = units.day(numpy.arange(20, 120 * 365.25, 12))
    timerange = Vq.arange(20 | units.day, 120 |units.yr, 10 |units.day)
    pdb.set_trace()

    instance = MercuryWayWard()
    instance.initialize_code()
    instance.central_particle.add_particles(sun)
    instance.orbiters.add_particles(planets)
    instance.commit_particles()

    channels = instance.orbiters.new_channel_to(planets)

    err = instance.evolve_model(10|units.day)
    pdb.set_trace()
    channels.copy()
    planets.savepoint(10|units.day)
    pdb.set_trace()

    for time in timerange:
        err = instance.evolve_model(time)
        channels.copy()
        planets.savepoint(time)

    instance.stop()

    pdb.set_trace()

if __name__ == "__main__":
    integrate_and_store()
