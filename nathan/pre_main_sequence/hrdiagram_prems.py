"""
Show how to generate Hertzsprung-Russell diagram
"""

import sys
import numpy
from matplotlib import pyplot
from amuse.plot import loglog, xlabel, ylabel

from amuse.units import units
from amuse.community.mesa.interface import MESA

from amuse import datamodel
from amuse.rfi.core import is_mpd_running
end_time = 2 | units.Gyr
stellar_mass =  2.0 | units.MSun

def simulate_evolution_tracks():
    stellar_evolution = MESA()#redirection="none")
    stellar_evolution.parameters.AGB_wind_scheme = 0
    stellar_evolution.parameters.RGB_wind_scheme = 0
#~    stellar_evolution.commit_parameters()

    stars = datamodel.Particles(6)
    stars.mass = [0.1, 0.5, 1.0, 3.0, 15.0, 45.0]|units.MSun
    stars = datamodel.Particles(3)
    stars.mass = [.25, .25, .25]|units.MSun
    stars = datamodel.Particles(1)
    stars.mass = 10.|units.MSun

    stars = stellar_evolution.pre_ms_stars.add_particles(stars)
#~    stars.mass_change = [1.0e-8, 1.0e-7, 1.0e-4] | units.MSun / units.yr
    stellar_evolution.commit_particles()
    
    luminosity_at_time = []
    temperature_at_time = []
    time = []
    
    for star in stars:
        one_luminosity_at_time = [] | units.LSun
        one_temperature_at_time = [] | units.K
        one_time = [] | units.yr
        while star.stellar_type == 17 | units.stellar_type:
            one_luminosity_at_time.append(star.luminosity)
            one_temperature_at_time.append(star.temperature)
            one_time.append(star.age)
            star.evolve_one_step()
#~            if star.age > 1e5 | units.yr:
#~                star.mass_change =  0.0 | units.MSun / units.yr
#~            print ".",
            print star.stellar_type, star.age, star.mass, star.luminosity, star.radius
        luminosity_at_time.append(one_luminosity_at_time)
        temperature_at_time.append(one_temperature_at_time)
        time.append(one_time)
        
    stellar_evolution.stop()
    
    return temperature_at_time, luminosity_at_time, time
    
def plot_track(temperature_at_time, luminosity_at_time):
    pyplot.figure(figsize = (8, 6))
    pyplot.title('Hertzsprung-Russell diagram', fontsize=12)
    
    print len(temperature_at_time)
    print len(luminosity_at_time)
    for temperature, luminosity in zip(temperature_at_time, luminosity_at_time):
        print len(temperature)
        print len(luminosity)
        loglog(temperature, luminosity, marker="s")
    xlabel('Effective Temperature')
    ylabel('Luminosity')
    pyplot.xlim(pyplot.xlim()[::-1])
    pyplot.ylim(.1,1.e4)
    pyplot.show()
    

if __name__ in ('__main__', '__plot__'):        
    temperatures, luminosities, time = simulate_evolution_tracks()
    plot_track(temperatures, luminosities)
    
