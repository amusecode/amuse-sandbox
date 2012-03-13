"""
   Visualization for simple N-body integration.
   Reads particle set from file (nbody.hdf5) and prints subsequent frames.
"""
import sys
import numpy
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel, plot
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser
from time import sleep
from numpy.random import uniform, seed
from amuse.community.hop.interface import Hop

def find_clumps(particles, unit_converter):
    
    hop = Hop(unit_converter)
    hop.particles.add_particles(particles)
    hop.calculate_densities()
    hop.do_hop()
    
    result = [x.get_intersecting_subset_in(particles) for x in hop.groups()]
    
    hop.stop()
    
    return result

def plot_clumps(groups, total_mass):
    number_of_particles_in_group = []
    fraction_of_mass_in_group =  []

    for group in groups:
        number_of_particles_in_group.append(len(group))
        fraction = (group.mass.sum()/total_mass)
        fraction_of_mass_in_group.append(fraction)
    
    figure = pyplot.figure(figsize= (12,6))
    
    
    subplot = figure.add_subplot(1, 2, 1)
    
    colormap = pyplot.cm.Paired
    for index, group in enumerate(groups):
        color = colormap(1.0 * index / len(groups))
        subplot.scatter(
            group.x.value_in(units.parsec),
            group.y.value_in(units.parsec),
            s = group.mass.value_in(units.MSun),
            edgecolors = color,
            facecolors = color
        )
    
    subplot.set_xlim(0,1)
    subplot.set_ylim(0,1)
    subplot.set_xlabel('x (parsec)')
    subplot.set_ylabel('y (parsec)')
    
    subplot = figure.add_subplot(1, 2, 2)
        
    subplot.plot(
        number_of_particles_in_group,
        fraction_of_mass_in_group,
    )
    
    subplot.set_xscale('log')
    subplot.set_yscale('log')
    subplot.set_xlabel('N')
    subplot.set_ylabel('df/d(Log_10 N)')
    
    figure.savefig('x.png')
    pyplot.show()

def main(filename = "gravhydro.hdf5", lim=-1):
    storage = store.StoreHDF(filename,"r")
    stars = storage.load()
    x = [] | units.AU
    y = [] | units.AU
    Mtot = stars[0].mass.sum()
    Rmax = stars[0].position.length().amax()
    unit_converter = nbody_system.nbody_to_si(Mtot, Rmax)
#    Rmin = 10| units.AU
    x0 = stars[0].x[0]
    y0 = stars[0].y[0]
    for si in stars.history:
        time = si.get_timestamp()
        x.append(si[0].x-x0)
        y.append(si[0].y-y0)
#        clumps = find_clumps(si, unit_converter)
#        for index, clump in enumerate(clumps):
#            print len(clump), clump.mass.sum()/Mtot,
#            scatter(
#                clump.x.value_in(units.AU),
#                clump.y.value_in(units.AU))

#    plot_clumps(clumps, Mtot)
#    scatter(x0, y0, s=100)
    pyplot.title("Fig 2a from Heemsker & Savonije")
    plot(x, y)
    xlabel("X")
    ylabel("Y")
    pyplot.show()

def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "gravhydro.hdf5",
                      help="output filename [gravhydro.hdf5]")
    result.add_option("-l", dest="lim", type="float", default = -1,
                      help="boxsize")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)



