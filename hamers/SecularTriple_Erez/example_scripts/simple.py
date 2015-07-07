from amuse.community.seculartriple_Erez.interface import SecularTriple
from amuse.datamodel import Particles
from amuse.units import units,constants,quantities

from matplotlib import pyplot

import numpy

CV_SUCCESS=0
CV_ROOT_RETURN=2
CV_WARNING=99

def make_triple_system():
    particles = Particles(5)
    masses = [1.2,1.0,0.5] | units.MSun
    N=len(masses)
    for index,mass in enumerate(masses):
        star = particles[index]
        star.is_binary = False
        star.mass = mass
        star.envelope_mass = 0.2 | units.MSun
        star.radius = 100.0 | units.RSun
        star.envelope_radius = 0.1 | units.RSun
        star.luminosity = 1.0 | units.LSun
        star.apsidal_motion_constant = 0.1
        star.gyration_radius = 0.1
        star.stellar_type = 1 | units.stellar_type
        star.spin_angular_frequency = 1.0e6 | 1.0/units.yr

    particles[3].is_binary = True ### inner binary
    particles[3].child1 = particles[0] ### primary
    particles[3].child2 = particles[1] ### secondary
    particles[3].semimajor_axis = 1.0 | units.AU
    particles[3].eccentricity = 0.1            
    particles[3].inclination = 80.0*numpy.pi/180.0
    particles[3].argument_of_pericenter = 0.01
    particles[3].longitude_of_ascending_node = 0.01
    
    particles[4].is_binary = True ### outer binary
    particles[4].child1 = particles[2] ### tertiary
    particles[4].child2 = particles[3] ### inner binary
    particles[4].semimajor_axis = 100.0 | units.AU
    particles[4].eccentricity = 0.5
    particles[4].inclination = 0.0
    particles[4].argument_of_pericenter = 1.51
    particles[4].longitude_of_ascending_node = 1.41
    return particles
    
def evolve_triple_system(particles,end_time,output_time_step):
    code = SecularTriple(redirection="none")
    code.particles.add_particles(particles)
    
    code.parameters.verbose = False
    code.parameters.include_quadrupole_terms = True
    code.parameters.include_octupole_terms = True
    code.parameters.include_1PN_inner_terms = False
    code.parameters.include_1PN_outer_terms = False
    code.parameters.include_25PN_inner_terms = False
    code.parameters.include_25PN_outer_terms = False
    code.parameters.include_inner_tidal_terms = False
    code.parameters.include_outer_tidal_terms = False
    
    code.parameters.check_for_dynamical_stability = True
    code.parameters.check_for_inner_collision = True
    code.parameters.check_for_outer_collision = False
    code.parameters.check_for_inner_RLOF = False
    code.parameters.check_for_outer_RLOF = False
    
    binaries = particles[particles.is_binary]
    stars = particles - binaries
    
    channel_from_particles_to_code = particles.new_channel_to(code.particles)
    channel_from_code_to_particles = code.particles.new_channel_to(particles)

    ### quantities later used for plotting ###
    times_array = quantities.AdaptingVectorQuantity() 
    a_in_array = quantities.AdaptingVectorQuantity()
    e_in_array = []
    i_tot_array = []
    g_in_array = []

    time = 0.0 | units.Myr
    while (time < end_time):
        code.evolve_model(time)  
        
        flag = code.flag
        if flag == CV_SUCCESS:
            pass ### succesful integration
        elif flag ==  CV_ROOT_RETURN:
            print 'root found at t/Myr = ',code.model_time,'; root_finding_flag',code.root_finding_flag
            if code.dynamical_instability == True: print 'dynamical_instability'
            if code.inner_collision == True: print 'inner_collision'
            if code.outer_collision == True: print 'outer_collision'
            if code.RLOF_star1 == True: print 'RLOF_star1'
            if code.RLOF_star2 == True: print 'RLOF_star2'
            if code.RLOF_star3 == True: print 'RLOF_star3'
            
            break
        elif flag == CV_WARNING:
            print 'warning occurred during integration'
        elif flag <= 0:
            print 'unrecoverable error occurred during integration'
            break
            
        time += output_time_step

        times_array.append(time)
        
        channel_from_code_to_particles.copy()
        e_in_array.append(binaries[0].eccentricity)
        a_in_array.append(binaries[0].semimajor_axis)
        g_in_array.append(binaries[0].argument_of_pericenter)    
        i_tot_array.append(binaries[0].inclination)

        print 'time/Myr',time.value_in(units.Myr),'ein',binaries[0].eccentricity,'ain/AU',binaries[0].semimajor_axis.value_in(units.AU),'Omega1/(/yr)',stars[0].spin_angular_frequency.value_in(1.0/units.yr)
    
       # print 'Edot1',stars[0].tidal_E_dot
        
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
    particles = make_triple_system()
    
    output_time_step = 0.01 | units.Myr
    end_time = 2.0 | units.Myr
    times_array,a_in_array,e_in_array,g_in_array,i_tot_array = evolve_triple_system(particles,end_time,output_time_step)
    plot_function(times_array,a_in_array,e_in_array,g_in_array)

