import sys, os
import numpy
import time as pytime
import platform
from numpy import random
from amuse.lab import *
from amuse.couple import bridge
from make_pebble_disk import make_pebble_disk
from amuse.units.optparse import OptionParser
from evolve_disk import *
from make_planets import *
from amuse.units import quantities

def scatter_plot(stars, time, nstep, subplot, lim = 50):
    MEarth = 5.97219e+24 * units.kg
    
    m = 100*stars.mass/stars.mass.max()
    colors = ['r', 'b', 'g']
    particles = ParticlesSuperset([stars, stars[0].planets, stars[0].disk_particles])
    
    particles.move_to_center()
    star = stars[0]
    offset = star.position
    stars.position -= offset
    
    subplot.scatter(stars.x.value_in(units.AU), stars.y.value_in(units.AU), c='r', s=m)
    
    for star, color in zip(stars, colors):
        if len(star.disk_particles) > 0:
            star.disk_particles.position -= offset
            subplot.scatter(star.disk_particles.x.value_in(units.AU), star.disk_particles.y.value_in(units.AU), c=color, s=1)
    
    colors = ['r', 'b', 'g','y']
    for star, color in zip(stars, colors):
        star.planets.position -= offset
        subplot.scatter(star.planets.x.value_in(units.AU), star.planets.y.value_in(units.AU), color=colors, s=10)

    subplot.set_xlabel("X [AU]")
    subplot.set_ylabel("Y [AU]")
    subplot.set_title("Time = "+str(time.as_quantity_in(units.yr))+ " I: " +str(nstep))

    subplot.set_xlim(-lim, lim)
    subplot.set_ylim(-lim, lim)
    

def ae_plot(particles, time, nstep, subplot, xlim = (0,20)):
    subplot.set_ylim(0, 1)
    subplot.semilogx()
    mscale = 0.1*particles.mass.min()
    if mscale < 10 | units.kg:
        scale = 10
    else:
        scale = particles.mass / mscale
    subplot.scatter(particles.semimajor_axis.value_in(units.AU), particles.eccentricity, c='r',  s=scale, lw=0)
    subplot.set_xlabel("a [AU]")
    subplot.set_ylabel("e")

    subplot.set_xlim(*xlim)
    subplot.set_ylim(0,1)
    
    
def get_kepler_elements(star, pebble, kepler_code):
    pos = star.position - pebble.position
    vel = star.velocity - pebble.velocity
    kepler_code.initialize_from_dyn(star.mass + pebble.mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2])
    return kepler_code.get_elements()
    

def determine_orbital_parameters_of_pebbles(stars, kepler_code):
    for star in stars:
        for pebble in star.disk_particles:
            if pebble.star is None:
                continue
            a, e = get_kepler_elements(star, pebble, kepler_code)
            pebble.semimajor_axis = a
            pebble.eccentricity = e
            
def determine_orbital_parameters_of_planets(stars, kepler_code):
    for star in stars:
        for planet in star.planets:
            if planet.star is None:
                continue
            a, e = get_kepler_elements(star, planet, kepler_code)
            planet.semimajor_axis = a
            planet.eccentricity = e


def new_option_parser():
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", default = "nice_model.h5",
                      help="input filename")
    result.add_option("--no-gui",
                  action="store_false", dest="with_gui", default=True,
                  help="don't show gui, just create png files")
    result.add_option("-t", "--plot-type", 
                      dest="plot_type", default = "xy",
                      help="plot type (xy, xz, planets, pebbles)")
    result.add_option("-i", "--start-index", 
                      dest="start_index", default = 0, type='int',
                      help="first index of the history to select  (default: start of history)")
    result.add_option("-j", "--end-index", 
                      dest="end_index", default = -1, type='int',
                      help="last index of the history to select (default: end of history)")
    return result


if __name__ == "__main__":
    o, arguments  = new_option_parser().parse_args()
    #    random.seed(seed=o.seed)
    
    if not o.filename or not os.path.exists(o.filename):
        print "Need to specifiy the name of the file to plot"
        sys.exit(1)
    stored_stars = read_set_from_file(o.filename, "amuse")
    history = list(stored_stars.iter_history())
    all_steps = list(reversed(history[1:]))
    converter = nbody_system.nbody_to_si(stored_stars.mass.sum(), 100 | units.AU)
    if o.end_index < 0:
        o.end_index = None
    all_steps = all_steps[o.start_index:o.end_index]
    
    time = 1 | units.yr
    
    if not o.with_gui:
        import matplotlib
        matplotlib.use('Agg')

    kepler_code = None
    if o.plot_type == 'planets' or o.plot_type == 'pebbles':
        kepler_code = Kepler(converter)
        kepler_code.initialize_code()
        
    from matplotlib import pyplot
    
    figure = pyplot.figure(figsize=(24,24))
    step_size = max(len(all_steps) / 9, 1)
    print len(all_steps) / step_size, step_size, len(all_steps)
    while len(all_steps) / step_size < 9:
        if step_size == 1:
            break
        step_size -= 1
    nstep = 1
    for i, x in enumerate(all_steps[::step_size]):
        nstep = o.start_index + (i * step_size)
        if (i > 8): 
            break
        stars = x.copy()
        if i == 0:
            print x.collection_attributes
        else:
            print x.collection_attributes.model_time.as_quantity_in(units.Myr)
        time = x.collection_attributes.model_time
        
        if len(stars[0].disk_particles) > 0:
            stars[0].disk_particles.star = stars[0]
        stars[0].planets.star = stars[0]
        subplot = figure.add_subplot(3,3,i+1)
        if o.plot_type == 'planets':
            determine_orbital_parameters_of_planets(stars, kepler_code)
            ae_plot(stars[0].planets, time, nstep, subplot)
        elif o.plot_type == 'pebbles':
            determine_orbital_parameters_of_pebbles(stars, kepler_code)
            ae_plot(stars[0].disk_particles, time, nstep, subplot, (10,100))
        else:
            scatter_plot(stars, time, nstep, subplot)
        nstep += 1
    
    if o.with_gui:
        pyplot.show()
        
    figure.savefig(o.filename + '-'+ o.plot_type +'.png')#,bbox_inches='tight')
#   mencoder "mf://xy*.png" -mf fps=20 -ovc x264 -o movie.avi

