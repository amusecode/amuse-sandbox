import os
import os.path
import shutil
import math
import subprocess
import numpy

from amuse.units import units, constants, nbody_system
from amuse.datamodel import Particles, Particle, ParticlesSuperset
from amuse.support.exceptions import AmuseException
from amuse.io import write_set_to_file, read_set_from_file
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH

from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.community.huayno.interface import Huayno
from amuse.community.mi6.interface import MI6
from amuse.community.mesa.interface import MESA as stellar_evolution_code

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
from amuse.plot import scatter, xlabel, ylabel, plot,loglog,semilogx,semilogy, sph_particles_plot
from amuse.plot import pynbody_column_density_plot, HAS_PYNBODY


def new_hydro(sph_code, sph_particles, core, t_end, n_steps, core_radius):
    unit_converter = nbody_system.nbody_to_si(sph_particles.total_mass() + core.mass, t_end)
    system = sph_code(unit_converter, redirection="file", redirect_file="sph_code_out.log")
    if sph_code is Gadget2:
        system.parameters.epsilon_squared = core_radius**2
        system.parameters.max_size_timestep = t_end / n_steps
        system.parameters.time_max = 1.1 * t_end
        system.parameters.time_limit_cpu = 70.0 | units.day
    else:
        system.parameters.timestep = t_end / n_steps
        system.parameters.eps_is_h_flag = True
        core.radius = core_radius * 2
    
    system.dm_particles.add_particle(core)
    system.gas_particles.add_particles(sph_particles)
    return system

def evolve_system(system, t_end, n_steps):
    times = (t_end * range(1, n_steps+1) / n_steps).as_quantity_in(units.day)
    
    for i_step, time in enumerate(times):
        system.evolve_model(time)
        print "   Evolved to:", time,
        
        figname = os.path.join("plots", "hydro_giant_{0:=04}.png".format(i_step))
        print "  -   Hydroplot saved to: ", figname
        pynbody_column_density_plot(system.gas_particles, width=30|units.AU, vmin=26, vmax=33)
        scatter(system.dm_particles.x, system.dm_particles.y, c="w")
        pyplot.savefig(figname)
        pyplot.close()
        
        file_base_name = os.path.join("snapshots", "hydro_giant_{0:=04}_".format(i_step))
        write_set_to_file(system.gas_particles, file_base_name+"gas.amuse", format='amuse')
        write_set_to_file(system.dm_particles, file_base_name+"dm.amuse", format='amuse')
    
    system.stop()

def new_working_directory():
    i = 0
    current_directory = os.getcwd()
    while os.path.exists(os.path.join(current_directory, "run_{0:=03}".format(i))):
        i += 1
    new_directory = os.path.join(current_directory, "run_{0:=03}".format(i))
    os.mkdir(new_directory)
    print "Created new directory for output:", new_directory
    os.mkdir(os.path.join(new_directory, "plots"))
    os.mkdir(os.path.join(new_directory, "snapshots"))
    shutil.copy(__file__, new_directory)
    os.chdir(new_directory)

if __name__ == "__main__":
    #~sph_code = Fi
    sph_code = Gadget2
    
    run_name = "200k_837RSun"
    
    file_base_name = os.path.join("initial_conditions", run_name, "hydro_triple_")
    gas = read_set_from_file(file_base_name+"gas.amuse", format='amuse')
    core = read_set_from_file(file_base_name+"core.amuse", format='amuse')[0]
    new_working_directory()
    
    t_end = 100.0 | units.day
    n_steps = 1000
    print "\nSetting up {0} to simulate triple system".format(sph_code.__name__)
    hydro = new_hydro(sph_code, gas, core, t_end, n_steps, core.radius)
    
    print "\nEvolving to {0}".format(t_end)
    evolve_system(hydro, t_end, n_steps)
    
    print "Done"
