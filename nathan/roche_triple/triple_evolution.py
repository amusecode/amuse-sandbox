import os
import os.path
import shutil
import math
import subprocess
import numpy

from amuse.units import units, constants
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.datamodel import Particles, Particle, ParticlesSuperset
from amuse.io import write_set_to_file

from amuse.ext.star_to_sph import convert_stellar_model_to_SPH
from amuse.community.evtwin.interface import EVtwin
from amuse.community.gadget2.interface import Gadget2

try:
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
    from amuse.plot import scatter, xlabel, ylabel, plot,loglog,semilogx,semilogy, sph_particles_plot
except ImportError:
    HAS_MATPLOTLIB = False

def inject_energy(gas_particles):
    print "Adding (specific internal) energy to each of them."
    gas_particles.u += 1E-5*gas_particles.potential_energy()/ gas_particles.total_mass() 

def Rate_Mass_Transfer(donor_mass,radius,roche_radius):
    orbital_period=144.543977075 | units.day
    transfer_rate=[] | (units.MSun/units.day)
    for i in range(len(donor_mass)):
        d_m=donor_mass[i]
        r=radius[i]
        r_l=roche_radius[i]
        delta_R=r-r_l
        mass_transfer=-(d_m/orbital_period)*((delta_R/r_l)**3)
        transfer_rate.append(mass_transfer)
    return transfer_rate




def new_working_directory():
    i = 0
    current_directory = os.getcwd()
    while os.path.exists(os.path.join(current_directory, "triple_run_{0:=03}".format(i))):
        i += 1
    new_directory = os.path.join(current_directory, "triple_run_{0:=03}".format(i))
    os.mkdir(new_directory)
    print "Created new directory for output:", new_directory
    os.mkdir(os.path.join(new_directory, "plots"))
    os.mkdir(os.path.join(new_directory, "snapshots"))
    shutil.copy(__file__, new_directory)
    os.chdir(new_directory)

def get_relative_velocity(total_mass, semimajor_axis, ecc):
    return (constants.G * total_mass * ((1.0 + ecc)/(1.0 - ecc)) / semimajor_axis).sqrt()

def set_up_initial_conditions():
    stars = set_up_inner_binary()
    giant = set_up_outer_star(stars.total_mass())
    view_on_giant = stars.add_particle(giant)
    stars.move_to_center()
    return stars, view_on_giant

def set_up_inner_binary():
    semimajor_axis = 0.133256133158 | units.AU
    eccentricity = 0
    masses = [3.2, 3.1] | units.MSun
    orbital_period = (4 * numpy.pi**2 * semimajor_axis**3 / 
        (constants.G * masses.sum())).sqrt().as_quantity_in(units.day)
    
    print "   Initializing inner binary"
    print "   Orbital period inner binary:", orbital_period
    stars =  Particles(2)
    stars.mass = masses
    stars.position = [0.0, 0.0, 0.0] | units.AU
    stars.velocity = [0.0, 0.0, 0.0] | units.km / units.s
    stars[0].x = semimajor_axis
    stars[0].vy = get_relative_velocity(stars.total_mass(), semimajor_axis, eccentricity)
    stars.move_to_center()
    return stars

def set_up_outer_star(inner_binary_mass):
    semimajor_axis = 1.22726535008 | units.AU
    eccentricity = 0.15
    inclination = math.radians(9.0)
    
    print "   Initializing outer star"
    giant = Particle()
    giant.mass = 5.5 | units.MSun
    giant.position = semimajor_axis * ([math.cos(inclination), 0, math.sin(inclination)] | units.none)
    giant.velocity = get_relative_velocity(giant.mass + inner_binary_mass, 
        semimajor_axis, eccentricity) * ([0, 1, 0] | units.none)
    return giant

def estimate_roche_radius(triple, view_on_giant):
    # 'mass ratio' of giant to inner binary
    q = (view_on_giant.mass / (triple-view_on_giant).total_mass()).value_in(units.none)
    # Assume ~ circular orbit:
    a = (view_on_giant.position - (triple-view_on_giant).center_of_mass()).length()
    q13 = q**(1./3.)
    q23 = q13**2
    return (a*(0.49*q23/(0.6*q23+math.log(1+q13)))).as_quantity_in(units.RSun)

def evolve_stars(triple, view_on_giant, stellar_evolution_code):
    stop_radius = estimate_roche_radius(triple, view_on_giant)
    stellar_evolution = stellar_evolution_code(redirection='file', redirect_file='stellar_evolution_code_out.log')
    se_stars = stellar_evolution.particles.add_particles(triple)
    view_on_se_giant = se_stars - (triple - view_on_giant)
    time = 0 | units.Myr
    while (view_on_se_giant.radius < stop_radius):
        time += 1 | units.Myr
        stellar_evolution.evolve_model(time)
#~        break
    
    return se_stars, stellar_evolution

def convert_stars_to_sph(se_stars, view_on_giant, gd_stars, giant_only = False):
    number_of_sph_particles = 15000 # total; if not giant_only: n(giant) + n(inner binary)
    
    view_on_se_binary = se_stars - view_on_giant
    view_on_se_giant = (se_stars - view_on_se_binary)[0]
    
    # Make sure that the sph particles are of (almost) equal mass
    if giant_only:
        n_sph_giant = number_of_sph_particles
    else:
        n_sph_giant = int(round(number_of_sph_particles * 
            (view_on_se_giant.mass / se_stars.total_mass()).value_in(units.none)))
        n_sph_ms1 = int(round((number_of_sph_particles - n_sph_giant) * 
            (view_on_se_binary[0].mass / view_on_se_binary.total_mass()).value_in(units.none)))
        n_sph_ms2 = number_of_sph_particles - n_sph_giant - n_sph_ms1
    
    print "   Converting stellar model of giant to", n_sph_giant, "SPH particles"
    giant_in_sph = convert_stellar_model_to_SPH(
        view_on_se_giant,
        n_sph_giant,
        seed = 12345, # Do I care?
        with_core_particle = True,
        target_core_mass  = 1.4 | units.MSun # Do I care?
    )
    print "   Giant core radius:", giant_in_sph.core_radius.as_quantity_in(units.RSun)

    sph_particles = Particles()
    giant_in_sph.gas_particles.position += view_on_giant.position
    giant_in_sph.gas_particles.velocity += view_on_giant.velocity
    sph_particles.add_particles(giant_in_sph.gas_particles)
    
    gd_particles = Particles()
    giant_in_sph.core_particle.position += view_on_giant.position
    giant_in_sph.core_particle.velocity += view_on_giant.velocity
    gd_particles.add_particle(giant_in_sph.core_particle)
    
    star_views = [ParticlesSuperset([giant_in_sph.gas_particles, giant_in_sph.core_particle.as_set()])]
    view_on_gd_binary = gd_stars - view_on_giant
    if giant_only:
        print "   Adding inner binary as two point particles"
        gd_particles.add_particles(view_on_gd_binary)
        star_views.extend([view_on_gd_binary[:1], view_on_gd_binary[1:]])
    else:
        print "   Converting stellar model of star in inner binary, with mass", 
        print view_on_se_binary[0].mass, " to", n_sph_ms1, "SPH particles"
        ms1_in_sph = convert_stellar_model_to_SPH(view_on_se_binary[0], n_sph_ms1)
        ms1_in_sph.gas_particles.position += view_on_gd_binary[0].position
        ms1_in_sph.gas_particles.velocity += view_on_gd_binary[0].velocity
        sph_particles.add_particles(ms1_in_sph.gas_particles)
        print "   Converting stellar model of star in inner binary, with mass", 
        print view_on_se_binary[1].mass, " to", n_sph_ms2, "SPH particles"
        ms2_in_sph = convert_stellar_model_to_SPH(view_on_se_binary[1], n_sph_ms2)
        ms2_in_sph.gas_particles.position += view_on_gd_binary[1].position
        ms2_in_sph.gas_particles.velocity += view_on_gd_binary[1].velocity
        sph_particles.add_particles(ms2_in_sph.gas_particles)
        star_views.extend([ms1_in_sph.gas_particles, ms2_in_sph.gas_particles])
    
    return sph_particles, gd_particles, giant_in_sph.core_radius, star_views
    
def evolve_gas(sph_particles, gd_particles, sph_code, core_radius, star_views):
    t_end = 90.0 | units.day
    n_steps = 1000
    hydro_code_options = dict(number_of_workers=2, redirection='file', redirect_file='hydrodynamics_code_out.log')
    
    unit_converter = ConvertBetweenGenericAndSiUnits(
        1.0 | units.RSun, 
        sph_particles.total_mass() + gd_particles.total_mass(), 
        t_end)
    
    hydrodynamics = sph_code(unit_converter, **hydro_code_options)
    hydrodynamics.parameters.epsilon_squared = core_radius**2
    hydrodynamics.parameters.max_size_timestep = t_end
    hydrodynamics.parameters.time_max = 1.1 * t_end
    hydrodynamics.parameters.time_limit_cpu = 7.0 | units.day
    hydrodynamics.gas_particles.add_particles(sph_particles)
    hydrodynamics.dm_particles.add_particles(gd_particles)
    
    giant_center_of_mass = [] | units.RSun
    ms1_center_of_mass = [] | units.RSun
    ms2_center_of_mass = [] | units.RSun
    giant_center_of_mass_velocity = [] | units.km / units.s
    ms1_center_of_mass_velocity = [] | units.km / units.s
    ms2_center_of_mass_velocity = [] | units.km / units.s
    
    print "   Evolving for", t_end
    times = (t_end * range(1, n_steps+1) / n_steps).as_quantity_in(units.day)#[:10]
    try:
        for i_step, time in enumerate(times):
            hydrodynamics.evolve_model(time)
            print "   Evolved to:", time,
            giant_center_of_mass.append(star_views[0].center_of_mass())
            ms1_center_of_mass.append(star_views[1].center_of_mass())
            ms2_center_of_mass.append(star_views[2].center_of_mass())
            giant_center_of_mass_velocity.append(star_views[0].center_of_mass_velocity())
            ms1_center_of_mass_velocity.append(star_views[1].center_of_mass_velocity())
            ms2_center_of_mass_velocity.append(star_views[2].center_of_mass_velocity())
    #~        sph_star_radius.append(outer_sph_star.virial_radius())
    #~        roche_radius.append(Roche_Lobe.Roche_Lobe_radius(semi_major_axis,3.2,3.1,5.5))
    #~        r_sroche=Roche_Lobe.Roche_Lobe_radius(semi_major_axis,3.2,3.1,5.5)
    #~        outer=outer_sph_star.select(lambda x, y, z : x**2 + y**2 + z**2 >r_sroche**2, ["x", "y", "z"])
    #~        print 'mass_out',outer.total_mass()
    #~        m_out=core.mass+outer_sph_star.total_mass()-outer.total_mass()
    #~        donor_mass.append(m_out)
    #~        m_ins=primary.total_mass()+secondary.total_mass()+outer.total_mass()
    #~        m_in.append(m_ins)
            
            if not i_step % 10:
                snapshotfile = os.path.join("snapshots", "hydro_triple_{0:=04}_gas.amuse".format(i_step))
                write_set_to_file(hydrodynamics.gas_particles, snapshotfile, format='amuse')
                snapshotfile = os.path.join("snapshots", "hydro_triple_{0:=04}_dm.amuse".format(i_step))
                write_set_to_file(hydrodynamics.dm_particles, snapshotfile, format='amuse')
            
            pyplot.figure(figsize = [10, 10])
            sph_particles_plot(hydrodynamics.gas_particles, gd_particles=hydrodynamics.dm_particles, 
                view=[-1.2, 1.2, -1.2, 1.2] | units.AU)
            figname = os.path.join("plots", "hydro_triple{0:=04}.png".format(i_step))
            pyplot.savefig(figname)
            print "  -   Hydroplot was saved to: ", figname
            pyplot.close()
        
        hydrodynamics.stop()
    finally:
        print "   Creating movie from snapshots"
        subprocess.call(['mencoder', "mf://hydro_triple*.png", '-ovc', 'lavc', 
            '-o', '../hydro_triple.avi', '-msglevel', 'all=1'], cwd="./plots")
        
        print "   Calculating semimajor axis and eccentricity evolution for inner binary"
        # Some temporary variables to calculate semimajor_axis and eccentricity evolution
        ms1_total_mass = star_views[1].total_mass()
        ms2_total_mass = star_views[2].total_mass()
        total_mass = ms1_total_mass + ms2_total_mass
        rel_position = ms1_center_of_mass - ms2_center_of_mass
        rel_velocity = ms1_center_of_mass_velocity - ms2_center_of_mass_velocity
        separation = rel_position.lengths()
        speed_squared = rel_velocity.lengths_squared()
        
        # Now calculate the important quantities:
        semimajor_axis_binary = (constants.G * total_mass * separation / 
            (2 * constants.G * total_mass - separation * speed_squared)).as_quantity_in(units.AU)
        eccentricity_binary = (1.0 - (rel_position.cross(rel_velocity)**2).sum(axis=1) / 
            (constants.G * total_mass * semimajor_axis_binary)).sqrt().value_in(units.none)
        
        print "   Calculating semimajor axis and eccentricity evolution of the giant's orbit"
        # Some temporary variables to calculate semimajor_axis and eccentricity evolution
        rel_position = (ms1_total_mass * ms1_center_of_mass + 
            ms2_total_mass * ms2_center_of_mass)/total_mass - giant_center_of_mass
        rel_velocity = (ms1_total_mass * ms1_center_of_mass_velocity + 
            ms2_total_mass * ms2_center_of_mass_velocity)/total_mass - giant_center_of_mass_velocity
        total_mass += star_views[0].total_mass()
        separation = rel_position.lengths()
        speed_squared = rel_velocity.lengths_squared()
        
        # Now calculate the important quantities:
        semimajor_axis_giant = (constants.G * total_mass * separation / 
            (2 * constants.G * total_mass - separation * speed_squared)).as_quantity_in(units.AU)
        eccentricity_giant = (1.0 - (rel_position.cross(rel_velocity)**2).sum(axis=1) / 
            (constants.G * total_mass * semimajor_axis_giant)).sqrt().value_in(units.none)
        
    #~    mass_transfer_rate= Rate_Mass_Transfer(donor_mass,sph_star_radius,roche_radius)    
    #~    plot_mass_transfer(mass_transfer_rate, star_age)
        orbit_parameters_plot(semimajor_axis_binary, semimajor_axis_giant, times[:len(semimajor_axis_binary)])
        orbit_ecc_plot(eccentricity_binary, eccentricity_giant, times[:len(eccentricity_binary)])


def plot_radius(radius,time,roche_radius):
    figure = pyplot.figure(figsize = (6, 10), dpi = 100)
    subplot = figure.add_subplot(1, 1, 1)
    plot(time.as_quantity_in(units.Myr),radius )
    plot(time.as_quantity_in(units.Myr),roche_radius ) 
    pyplot.minorticks_on()
    xlabel('t')
    ylabel('Radius') 
    pyplot.text(20, 230, "Roche Lobe")
    pyplot.minorticks_on()
    pyplot.savefig('triple_outer_radius.png')

def plot_mass_transfer(mass_transfer_rate, times):
    figure = pyplot.figure(figsize = (8, 10), dpi = 100)
    subplot = figure.add_subplot(1, 1, 1)
    plot(times,mass_transfer_rate.as_quantity_in(units.MSun/units.day) )
    pyplot.minorticks_on()
    xlabel('t')
    ylabel('dM/dt') 
    pyplot.minorticks_on()
    pyplot.savefig('triple_mass_transfer.png')


def orbit_ecc_plot(eccentricity_in,eccentricity_out,time):
    figure = pyplot.figure(figsize = (10, 6), dpi = 100)
    subplot = figure.add_subplot(2, 1, 1)
    plot(time,eccentricity_in)
    xlabel('t')
    ylabel('e$_\mathrm{binary}$')
    
    subplot = figure.add_subplot(2, 1, 2)
    plot(time,eccentricity_out )  
    xlabel('t')
    ylabel('e$_\mathrm{giant}$')
    pyplot.minorticks_on()
    pyplot.savefig("triple_eccentricities.png")


def orbit_parameters_plot(semi_major_in,semi_major_out, time):
    figure = pyplot.figure(figsize = (10, 6), dpi = 100)
    subplot = figure.add_subplot(2, 1, 1)
    plot(time,semi_major_in )
    xlabel('t')
    ylabel('$a_\mathrm{binary}$')
    
    subplot = figure.add_subplot(2, 1, 2)
    plot(time,semi_major_out )
    xlabel('t')
    ylabel('$a_\mathrm{giant}$')
    pyplot.minorticks_on()
    pyplot.savefig("triple_semimajor_axis.png")


if __name__ == "__main__":
    stellar_evolution_code = EVtwin
    sph_code = Gadget2
    
    new_working_directory()
    
    print "Initializing triple"
    triple, view_on_giant = set_up_initial_conditions()
    print "\nInitialization done:\n", triple
    
    print "\nEvolving with", stellar_evolution_code.__name__
    se_stars, se_code_instance = evolve_stars(triple, view_on_giant, stellar_evolution_code)
    print "\nStellar evolution done:\n", se_stars
    
    print "\nConverting evolved stars to SPH particles"
    sph_particles, gd_particles, core_radius, star_views = convert_stars_to_sph(
        se_stars, 
        view_on_giant, 
        triple, 
        giant_only = False
    )
    se_code_instance.stop()
    
    print "Evolving with", sph_code.__name__
    evolve_gas(sph_particles, gd_particles, sph_code, core_radius, star_views)
    print "Done"
