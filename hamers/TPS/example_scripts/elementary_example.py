from amuse.community.seculartriple.interface import SecularTriple
from amuse.datamodel import Particles
from amuse.units import units,constants,quantities

from matplotlib import pyplot

import numpy

def make_triple_system():
    binaries = Particles(2)
    inner_stars = Particles(2)
    outer_stars = Particles(2)

    binaries[0].child1 = inner_stars[0]
    binaries[0].child2 = inner_stars[1]
    binaries[1].child1 = outer_stars[0]
    binaries[1].child2 = outer_stars[1]

    binaries[0].child1.mass = 1.0 | units.MSun
    binaries[0].child2.mass = 0.5 | units.MSun
    binaries[1].child1.mass = 0.5 | units.MSun
    binaries[0].child1.radius = 1.0 | units.RSun
    binaries[0].child2.radius = 1.0 | units.RSun
    binaries[1].child1.radius = 1.0 | units.RSun

    binaries[0].semimajor_axis = 1.0 | units.AU
    binaries[1].semimajor_axis = 100.0 | units.AU
    binaries[0].eccentricity = 0.1
    binaries[1].eccentricity = 0.5
    binaries[0].argument_of_pericenter = 0.0
    binaries[1].argument_of_pericenter = 0.0
    binaries[0].inclination = 80.0*numpy.pi/180.0 ### this is the inclination between the orbital planes of binaries[0] and binaries[1] ###

    return binaries
    
def evolve_triple_system(binaries,end_time,output_time_step):
    code = SecularTriple()
    code.binaries.add_particles(binaries)

    ### quantities later used for plotting ###
    times_array = quantities.AdaptingVectorQuantity() 
    e_in_array = []
    i_tot_array = []
    g_in_array = []

    time = 0.0 | units.Myr
    while (time < end_time):
        code.evolve_model(time)  
        time += output_time_step

        times_array.append(time)
        e_in_array.append(code.binaries[0].eccentricity)
        g_in_array.append(code.binaries[0].argument_of_pericenter)    
        i_tot_array.append(code.binaries[0].inclination)

        print 'time/Myr',time.value_in(units.Myr),'ein',code.binaries[0].eccentricity
    
    e_in_array = numpy.array(e_in_array)
    g_in_array = numpy.array(g_in_array)
    i_tot_array = numpy.array(i_tot_array)
    return times_array,e_in_array,g_in_array,i_tot_array
    
def plot_function(times_array,e_in_array,g_in_array):

    figure = pyplot.figure(figsize=(10,13))
    N_subplots = 3

    plot_e_in = figure.add_subplot(N_subplots,1,1)
    plot_i_tot = figure.add_subplot(N_subplots,1,2)
    plot_e_in_g_in = figure.add_subplot(N_subplots,1,3)

    times_array_Myr = times_array.value_in(units.Myr)
    plot_e_in.plot(times_array_Myr,e_in_array)
    plot_i_tot.plot(times_array_Myr,i_tot_array*180.0/numpy.pi)
    plot_e_in_g_in.plot(numpy.cos(g_in_array),e_in_array)

    t_max_Myr = max(times_array_Myr)
    plot_e_in.set_xlim(0,t_max_Myr)
    plot_i_tot.set_xlim(0,t_max_Myr)

    plot_i_tot.set_ylim(0.9*min(i_tot_array*180.0/numpy.pi),1.1*max(i_tot_array*180.0/numpy.pi))

    plot_e_in.set_xlabel('$t/\mathrm{Myr}$')
    plot_i_tot.set_xlabel('$t/\mathrm{Myr}$')
    plot_e_in_g_in.set_xlabel('$\cos(g_\mathrm{in})$')

    plot_e_in.set_ylabel('$e_\mathrm{in}$')
    plot_i_tot.set_ylabel('$i_\mathrm{tot} ({}^\circ)$')
    plot_e_in_g_in.set_ylabel('$e_\mathrm{in}$')
    figure.subplots_adjust(left=0.2, right=0.85, top=0.8, bottom=0.15)

    pyplot.show()

if __name__ in ('__main__','__plot__'):
    binaries = make_triple_system()
    
    output_time_step = 0.01 | units.Myr
    end_time = 10.0 | units.Myr
    times_array,e_in_array,g_in_array,i_tot_array = evolve_triple_system(binaries,end_time,output_time_step)
    plot_function(times_array,e_in_array,g_in_array)

