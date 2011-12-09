"""
AMUSE Tutorial, step A: equal-mass n-body dynamics
"""
import os
import os.path
import shutil
import subprocess
import pickle

from amuse.units import units, constants, nbody_system
from amuse.io import write_set_to_file, read_set_from_file
from amuse.ic.plummer import new_plummer_model

from amuse.community.hermite0.interface import Hermite

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel, plot

def new_working_directory():
    i = 0
    current_directory = os.getcwd()
    while os.path.exists(os.path.join(current_directory, "A_run_{0:=03}".format(i))):
        i += 1
    new_directory = os.path.join(current_directory, "A_run_{0:=03}".format(i))
    os.mkdir(new_directory)
    print "Created new directory for output:", new_directory
    os.mkdir(os.path.join(new_directory, "plots"))
    os.mkdir(os.path.join(new_directory, "snapshots"))
    shutil.copy(__file__, new_directory)
    os.chdir(new_directory)

def continue_evolution(dynamics_code, unit_converter, t_end, n_steps, number_of_workers):
    print "Loading snapshots...",
    files = os.listdir("snapshots")
    files.sort()
    files = files[-2:]
    print files
    stars = read_set_from_file(os.path.join("snapshots", files[0]), format='amuse')
    
    print "\nSetting up {0} system with a cluster of {1} stars".format(dynamics_code.__name__, len(stars))
    dynamics_system = set_up_dynamics_code(dynamics_code, unit_converter, stars, number_of_workers)
    
    print "\nEvolving for", t_end
    evolve_cluster(dynamics_system, t_end, n_steps, previous_data=os.path.join("snapshots", files[1]))
    print "Done"

def set_up_dynamics_code(dynamics_code, unit_converter, stars, number_of_workers):
    dynamics_system = dynamics_code(unit_converter, number_of_workers=number_of_workers, 
        redirection='file', redirect_file='dynamics_code_out.log')
    dynamics_system.parameters.epsilon_squared = 0 | units.m**2
    dynamics_system.particles.add_particles(stars)
    return dynamics_system

def evolve_cluster(dynamics_system, t_end, n_steps, previous_data=None):
    times = (t_end * range(1, n_steps+1) / n_steps).as_quantity_in(units.Myr)
    if previous_data:
        with open(previous_data, 'rb') as file:
            (all_times, potential_energies, kinetic_energies) = pickle.load(file)
        all_times.extend(times + all_times[-1])
    else:
        all_times = times
        potential_energies = dynamics_system.potential_energy.as_vector_with_length(1).as_quantity_in(units.erg)
        kinetic_energies = dynamics_system.kinetic_energy.as_vector_with_length(1).as_quantity_in(units.erg)
    i_offset = len(potential_energies) - 1
    
    for i_step, time in enumerate(times):
        dynamics_system.evolve_model(time)
        print "   Evolved to:", time,
        potential_energies.append(dynamics_system.potential_energy)
        kinetic_energies.append(dynamics_system.kinetic_energy)
        make_plot(dynamics_system.particles.position, i_step + i_offset, all_times[i_step + i_offset])
        
        if i_step % 10 == 9:
            snapshotfile = os.path.join("snapshots", "star_cluster_{0:=04}.amuse".format(i_step + i_offset))
            write_set_to_file(dynamics_system.particles, snapshotfile, format='amuse')
            datafile = os.path.join("snapshots", "star_cluster_{0:=04}_info.pkl".format(i_step + i_offset))
            with open(datafile, 'wb') as outfile:
                pickle.dump((all_times[:len(potential_energies)-1], potential_energies, 
                    kinetic_energies), outfile)
    
    dynamics_system.stop()
    
    print "   Creating movie from snapshots"
    subprocess.call(['mencoder', "mf://star_cluster_*.png", '-ovc', 'lavc', 
        '-o', '../star_cluster.avi', '-msglevel', 'all=1'], cwd="./plots")
    
    energy_evolution_plot(all_times, kinetic_energies, potential_energies)

def make_plot(positions, plot_number, time):
    x = positions.x.as_quantity_in(units.parsec)
    y = positions.y.as_quantity_in(units.parsec)
    pyplot.figure(figsize = [10, 10])
    pyplot.title("time = {0}".format(time))
    current_axes = pyplot.gca()
    current_axes.set_axis_bgcolor('#101010')
    current_axes.set_aspect("equal", adjustable = "box")
    current_axes.set_xlim(-2.0, 2.0, emit=True, auto=False)
    current_axes.set_ylim(-2.0, 2.0, emit=True, auto=False)
    scatter(x, y, c='w', marker='o')
    xlabel('x')
    ylabel('y')
    figname = os.path.join("plots", "star_cluster_{0:=04}.png".format(plot_number))
    pyplot.savefig(figname)
    print "  -   Plot was saved to: ", figname
    pyplot.close()

def energy_evolution_plot(time, kinetic, potential, figname = "energy_evolution.png"):
    time.prepend(0.0 | units.day)
    pyplot.figure(figsize = (5, 5))
    plot(time, kinetic, label='K')
    plot(time, potential, label='U')
    plot(time, kinetic + potential, label='E')
    xlabel('Time')
    ylabel('Energy')
    pyplot.legend(prop={'size':"x-small"}, loc=4)
    pyplot.savefig(figname)
    pyplot.close()

if __name__ == "__main__":
    dynamics_code = Hermite
    number_of_particles = 1000
    total_mass = 1000.0 | units.MSun
    virial_radius = 0.5 | units.parsec
    t_end = 10.0 | units.Myr
    n_steps = 1000
    number_of_workers = 4
    
    unit_converter = nbody_system.nbody_to_si(total_mass, virial_radius)
    
    if os.path.exists("snapshots"):
        print "Found snapshots folder, continuing evolution of previous run"
        continue_evolution(dynamics_code, unit_converter, t_end, n_steps, number_of_workers)
        exit(0)
    
    new_working_directory()
    
    print "\nGenerating initial conditions"
    stars = new_plummer_model(number_of_particles, unit_converter)
    
    print "\nSetting up {0} system with a cluster of {1} stars".format(dynamics_code.__name__, number_of_particles)
    dynamics_system = set_up_dynamics_code(dynamics_code, unit_converter, stars, number_of_workers)
    
    print "\nEvolving for {0} (n-body units: {1})".format(t_end, unit_converter.to_nbody(t_end))
    evolve_cluster(dynamics_system, t_end, n_steps)
    print "Done"
