import numpy
from amuse.community.sse.interface import SSE
from amuse.units import units
from amuse.ext.solarsystem import Solarsystem
from amuse.plot import *

from amuse import datamodel
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

class SSEWithMassEvolve(SSE):
    def __init__(self, **options):
        SSE.__init__(self, convert_nbody = None, **options)

    def evolve_mass(self, mass):
        timestep = 1.0 | units.Myr
        current_mass = self.particles[0].mass
        current_time = self.particles[0].age

        while current_mass >= mass:
            current_time += timestep
            self.evolve_model(current_time)
            current_mass = self.particles[0].mass
        return current_time

def plottillagb():
    sse = SSEWithMassEvolve()
    sse.commit_parameters()
    sun, planets = Solarsystem.new_solarsystem()
    sse.particles.add_particles(sun)
    sse.commit_particles()
    channel = sse.particles.new_channel_to(sun)
    channel.copy()
    
    massrange = units.MSun(numpy.arange(1, 0.8 ,-0.001))
    print massrange
    masses = []|units.MSun

    timerange = [] | units.Myr

    for mass in massrange:
        #sse.evolve_model(time)
        sse.evolve_mass(mass)
        timerange.append(sse.evolve_mass(mass))
        channel.copy()
        masses.append(sse.particles[0].mass)
        #print time.value_in(units.yr), sse.particles[0].mass.value_in(units.MSun)
        
    sse.stop()
    plot(massrange, timerange,'.')
    native_plot.show()

if __name__ == "__main__":
    plottillagb()
