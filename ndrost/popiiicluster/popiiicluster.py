import os.path
import numpy.random

from amuse.units import nbody_system, units
from amuse.io import write_set_to_file

from amuse.ic.kingmodel import new_king_model
from amuse.ic.flatimf import new_flat_mass_distribution

from amuse.couple.collision_handler import CollisionHandler
from amuse.couple.parallel_stellar_evolution import ParallelStellarEvolution
from amuse.support.project import new_working_directory

from amuse.community.ph4.interface import ph4
from amuse.community.mesa.interface import MESA
from amuse.community.mmams.interface import MakeMeAMassiveStar

def new_gravity_code(stars, convert_nbody):
    gravity = ph4(convert_nbody, mode='normal', 
        redirection='file', redirect_file=os.path.join('logs', 'ph4.log'))
    gravity.parameters.epsilon_squared = 1.0e-10 | nbody_system.length**2
    gravity.particles.add_particles(stars)
    gravity.commit_particles()
    gravity.stopping_conditions.collision_detection.enable()
    return gravity

def new_stellar_evolution_code(stars):
    number_of_workers = 2
    individual_options = [dict(
        redirect_file=os.path.join("logs", "MESA_{0}.log".format(i))
    ) for i in range(number_of_workers)]
    
    stellar_evolution = ParallelStellarEvolution(MESA, 
	must_run_threaded = False,
        number_of_workers=number_of_workers, 
        redirection='file', 
        individual_options=individual_options
    )
    stellar_evolution.parameters.metallicity = 0.0
    stellar_evolution.parameters.RGB_wind_scheme = 0
    stellar_evolution.parameters.AGB_wind_scheme = 0
    stellar_evolution.particles.add_particles(stars)
    return stellar_evolution

def new_collision_handler(stellar_evolution, gravity):
    collision_code = MakeMeAMassiveStar(redirection='file', redirect_file=os.path.join('logs', 'MakeMeAMassiveStar.log'))
    collision_code.parameters.target_n_shells_mixing = 2000
    collision_code.parameters.dump_mixed_flag = True
    collision_code.commit_parameters()
    handler = CollisionHandler(
        collision_code, 
        stellar_evolution_code = stellar_evolution,
        gravity_code = gravity,
        verbose = True
    )
    return handler

def new_stars(masses, convert_nbody):
    stars = new_king_model(len(masses), 6, convert_nbody)
    stars.mass = masses
    stars.scale_to_standard(convert_nbody=convert_nbody)
    return stars

def write_data(gravity, stellar_evolution, filename):
    write_set_to_file(
        gravity.particles, 
        filename, 'amuse', 
        gravity_time=gravity.model_time, 
        stellar_time=stellar_evolution.model_time)

def diagnostics(gravity, E0):
    print "   Evolved system (currently {0} particles) up to: {1}".format(len(gravity.particles), gravity.model_time.as_quantity_in(units.Myr))
    K = gravity.kinetic_energy
    U = gravity.potential_energy
    print "   Total energy E=K+U: {0}, ratio K/U: {1}, energy error dE/E0: {2}".format(K+U, K/U, (K+U-E0)/E0)

def safe_cleanup(*args):
    for code in args:
        try:
            code.stop()
        except:
            pass

def simulate_pop_III_cluster(masses, virial_radius, time_end, delta_t):
    convert_nbody = nbody_system.nbody_to_si(masses.sum(), virial_radius)
    stars = new_stars(masses, convert_nbody)
    gravity = new_gravity_code(stars, convert_nbody)
    collision_detection = gravity.stopping_conditions.collision_detection
    stellar_evolution = new_stellar_evolution_code(stars)
    collision_handler = new_collision_handler(stellar_evolution, gravity)
    
    from_stellar_evolution_to_gravity = stellar_evolution.particles.new_channel_to(gravity.particles)
    from_stellar_evolution_to_gravity.copy_attributes(["mass", "radius"])
    
    initial_total_energy = stars.kinetic_energy() + stars.potential_energy()
    output_file = 'popiiicluster_N{0}_history.hdf5'.format(len(masses))
    
    try:
        while gravity.model_time < time_end:
            print
            target_time = gravity.model_time + delta_t
            gravity_request = gravity.evolve_model.async(target_time)
            stellar_evolution.evolve_model(target_time)
            print "   Stellar evolution done"
            gravity_request.result()
            gravity.synchronize_model()
            print "   Gravity done, {0} collision(s) detected".format(len(collision_detection.particles(0)))
            from_stellar_evolution_to_gravity.copy_attributes(["mass", "radius"])
            write_data(gravity, stellar_evolution, output_file)
            if collision_detection.is_set():
                collision_handler.handle_collisions(collision_detection.particles(0), collision_detection.particles(1))
                write_data(gravity, stellar_evolution, output_file)
            diagnostics(gravity, initial_total_energy)
    finally:
        safe_cleanup(gravity, stellar_evolution, collision_handler.collision_code)
    

if __name__ == '__main__':
    new_working_directory(__file__, sub_directories=["logs"])
    
    number_of_stars = 46
    virial_radius = 0.1 | units.parsec
    
    time_end = 10.0 | units.Myr
    delta_t = 0.001 | units.Myr
    
    numpy.random.seed(54321)
    masses = new_flat_mass_distribution(number_of_stars, mass_min=1.0|units.MSun, mass_max=100|units.MSun).sorted()
    simulate_pop_III_cluster(masses, virial_radius, time_end, delta_t)
    print "Done!"
