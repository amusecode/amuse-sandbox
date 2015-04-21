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

    binaries[0].child1.mass = 1.3 | units.MSun
    binaries[0].child2.mass = 0.5 | units.MSun
    binaries[1].child1.mass = 0.5 | units.MSun
    binaries[0].child1.envelope_mass = 0.8 | units.MSun
    binaries[0].child2.envelope_mass = 0.01 | units.MSun
    binaries[1].child1.envelope_mass = 0.01 | units.MSun    
    binaries[0].child1.radius = 100.0 | units.RSun
    binaries[0].child2.radius = 1.0 | units.RSun
    binaries[1].child1.radius = 1.0 | units.RSun  
    binaries[0].child1.luminosity = 1.0 | units.LSun
    binaries[0].child2.luminosity = 1.0 | units.LSun
    binaries[1].child1.luminosity = 1.0 | units.LSun          
    binaries[1].child1.envelope_radius = 0.95 | units.RSun
    binaries[0].child1.envelope_radius = 0.01 | units.RSun
    binaries[0].child2.envelope_radius = 0.01 | units.RSun
    binaries[0].child1.apsidal_motion_constant = 0.1
    binaries[0].child2.apsidal_motion_constant = 0.1
    binaries[1].child1.apsidal_motion_constant = 0.1
    binaries[0].child1.gyration_radius = 0.1
    binaries[0].child2.gyration_radius = 0.1
    binaries[1].child1.gyration_radius = 0.1
    binaries[0].child1.stellar_type = 1
    binaries[0].child2.stellar_type = 1
    binaries[1].child1.stellar_type = 1
    binaries[0].child1.spin_angular_frequency = 1.0e3 | 1.0/units.yr
    binaries[0].child2.spin_angular_frequency = 1.0e3 | 1.0/units.yr
    binaries[1].child1.spin_angular_frequency = 1.0e3 | 1.0/units.yr       
    binaries[0].semimajor_axis = 1.0 | units.AU
    binaries[1].semimajor_axis = 100.0 | units.AU
    
    binaries[0].child1.mass_transfer_rate = -1.0e-6 | units.MSun / units.yr
    binaries[0].mass_transfer_accretion_parameter = 1.0
    binaries[1].mass_transfer_accretion_parameter = 1.0
    binaries[0].mass_transfer_angular_momentum_loss_parameter = 0.0
    binaries[1].mass_transfer_angular_momentum_loss_parameter = 0.0
    
    binaries[0].eccentricity = 0.1
    binaries[1].eccentricity = 0.5
    binaries[0].argument_of_pericenter = 0.5
    binaries[1].argument_of_pericenter = 0.1
    binaries[0].longitude_of_ascending_node = 0.1
    binaries[1].longitude_of_ascending_node = 0.3

    binaries[0].inclination = 80.0*numpy.pi/180.0 ### this is the inclination between the orbital planes of binaries[0] and binaries[1] ###
    binaries[1].inclination = 0.0

    return binaries
    
def evolve_triple_system(binaries,end_time,output_time_step):
    code = SecularTriple()
    code.binaries.add_particles(binaries)
    code.parameters.equations_of_motion_specification = 0
    code.parameters.f_quad = 1.0
    code.parameters.f_oct = 1.0    
    code.parameters.f_mass_transfer = 0.0
    code.parameters.f_1PN_in = 0.0
    code.parameters.f_1PN_out = 0.0    
    code.parameters.f_25PN_in = 0.0        
    code.parameters.f_25PN_out = 0.0        

#    code.parameters.f_quad = 1.0
#    code.parameters.f_oct = 1.0    
#    code.parameters.f_tides = 1.0
#    code.parameters.f_mass_transfer = 0.0

    ### quantities later used for plotting ###
    times_array = quantities.AdaptingVectorQuantity() 
    a_in_array = quantities.AdaptingVectorQuantity()
    e_in_array = []
    i_tot_array = []
    g_in_array = []

    time = 0.0 | units.Myr
    while (time < end_time):
        code.evolve_model(time)  
        time += output_time_step

        times_array.append(time)
        e_in_array.append(code.binaries[0].eccentricity)
        a_in_array.append(code.binaries[0].semimajor_axis)
        g_in_array.append(code.binaries[0].argument_of_pericenter)    
        i_tot_array.append(code.binaries[0].inclination)

        print 'time/Myr',time.value_in(units.Myr),'ein',code.binaries[0].eccentricity,'ain/AU',code.binaries[0].semimajor_axis.value_in(units.AU),'Omega1/(/yr)',code.binaries[0].child1.spin_angular_frequency.value_in(1.0/units.yr)
    
    e_in_array = numpy.array(e_in_array)
    g_in_array = numpy.array(g_in_array)
    i_tot_array = numpy.array(i_tot_array)
    return times_array,a_in_array,e_in_array,g_in_array,i_tot_array
    
def plot_function(times_array,a_in_array,e_in_array,g_in_array):

    figure = pyplot.figure(figsize=(10,13))
    N_subplots = 4

    plot_e_in = figure.add_subplot(N_subplots,1,1)
    plot_i_tot = figure.add_subplot(N_subplots,1,2)
    plot_e_in_g_in = figure.add_subplot(N_subplots,1,3)
    plot_a_in = figure.add_subplot(N_subplots,1,4)

    times_array_Myr = times_array.value_in(units.Myr)
    a_in_array_AU = a_in_array.value_in(units.AU)
    plot_e_in.plot(times_array_Myr,e_in_array)
    plot_i_tot.plot(times_array_Myr,i_tot_array*180.0/numpy.pi)
    plot_e_in_g_in.plot(numpy.cos(g_in_array),e_in_array)
    plot_a_in.plot(times_array_Myr,a_in_array_AU)

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
    times_array,a_in_array,e_in_array,g_in_array,i_tot_array = evolve_triple_system(binaries,end_time,output_time_step)
    plot_function(times_array,a_in_array,e_in_array,g_in_array)

