import os
import os.path
import shutil
import numpy

from amuse.units import nbody_system, units
from amuse.datamodel import Particle, ParticlesSuperset
from amuse.support.exceptions import AmuseException
from amuse.couple.parallel_stellar_evolution import ParallelStellarEvolution
from amuse.io import write_set_to_file

from amuse.community.ph4.interface import ph4
from amuse.community.hermite0.interface import Hermite
from amuse.community.mesa.interface import MESA
from amuse.community.sse.interface import SSE
from amuse.community.mmams.interface import MakeMeAMassiveStar

from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.flatimf import new_flat_mass_distribution


def new_working_directory():
    i = 0
    current_directory = os.getcwd()
    while os.path.exists(os.path.join(current_directory, "run_{0:=03}".format(i))):
        i += 1
    new_directory = os.path.join(current_directory, "run_{0:=03}".format(i))
    os.mkdir(new_directory)
    print "Created new directory for output:", new_directory
    os.mkdir(os.path.join(new_directory, "logs"))
    shutil.copy(__file__, new_directory)
    os.chdir(new_directory)

def new_gravity_code(stars, convert_nbody):
    gravity = ph4(convert_nbody, mode='gpu', channel_type='ibis', hostname='lgm19', 
        redirection='file', redirect_file=os.path.join('logs', 'gravity_code_out.log'))
#~    gravity.stopping_conditions.collision_detection.disable()
    gravity.stopping_conditions.collision_detection.enable()
    gravity.particles.add_particles(stars)
    return gravity

def new_stellar_evolution_codes(stars):
    n_low_mass = (stars.mass < 10|units.MSun).sum()
    print n_low_mass, "particles in SSE,", len(stars)-n_low_mass, "particles in MESA"
    
    number_of_workers = 54
    nodes = ['lgm00', 'lgm02', 'lgm03', 'lgm04', 'lgm05', 'lgm06', 'lgm07', 'lgm08', 'lgm09', 
        'lgm10', 'lgm11', 'lgm12', 'lgm13', 'lgm14', 'lgm15', 'lgm16', 'lgm17', 'lgm19']
    
    n_nodes = len(nodes)
    individual_options = [dict(
        hostname=nodes[i%n_nodes], 
        redirect_file=os.path.join("logs", "stellar_evolution_code_out_{0}.log".format(i))
    ) for i in range(number_of_workers)]
    
#~    stellar_evolution = ParallelStellarEvolution(MESA, number_of_workers=number_of_workers, 
#~        redirection='file', individual_options=individual_options)
    
    stellar_evolution = ParallelStellarEvolution(MESA, number_of_workers=number_of_workers, 
        channel_type='ibis', redirection='file', individual_options=individual_options,
        output_data_root_directory="/home/vriesn/amuse-svn/data", 
        input_data_root_directory="/home/vriesn/amuse-svn/data", 
        default_path_to_MESA_data="/home/vriesn/amuse-svn/src/amuse/community/mesa/src/mesa/data")
    
#~    stellar_evolution = MESA(channel_type='ibis', hostname='lgm11', redirection='file', redirect_file="stellar_evolution_code_out_0.log",
#~        output_data_root_directory="/home/vriesn/amuse-svn/data", 
#~        input_data_root_directory="/home/vriesn/amuse-svn/data", 
#~        default_path_to_MESA_data="/home/vriesn/amuse-svn/src/amuse/community/mesa/src/mesa/data")
    
#~    stellar_evolution = MESA(redirection='file', redirect_file="stellar_evolution_code_out_0.log")
    
    stellar_evolution.parameters.metallicity = 0.0
    stellar_evolution.particles.add_particles(stars[n_low_mass:])
    
    se_light = SSE(redirection='file', redirect_file=os.path.join("logs", "low_mass_stellar_evolution_code_out.log"))
    se_light.parameters.metallicity = 0.0001
    se_light.particles.add_particles(stars[:n_low_mass])
    return stellar_evolution, se_light

def move_from_light_to_real_stellar_evolution(missing, se_light, stellar_evolution):
    print "   Move the missing particles from the 'light' code to the 'real' stellar evolution code:\n", missing
    vivified = stellar_evolution.particles.add_particles(missing)
    for star in vivified:
        star.evolve_for(stellar_evolution.model_time)
    se_light.particles.remove_particles(missing)
    print vivified

def main():
    number_of_stars = 35000
    time_end = 10.0 | units.Myr
    delta_t = 0.0001 | units.Myr

    numpy.random.seed(12345)
    masses = new_salpeter_mass_distribution(number_of_stars, 
#~    masses = new_flat_mass_distribution(number_of_stars, 
        mass_min = 0.1 | units.MSun, mass_max = 100 | units.MSun).sorted()
    
    convert_nbody = nbody_system.nbody_to_si(masses.sum(), 1.0 | units.parsec)
    stars = new_plummer_model(number_of_stars, convert_nbody)
    stars.mass = masses
    stars.scale_to_standard(convert_nbody=convert_nbody)

    gravity = new_gravity_code(stars, convert_nbody)
    collision_detection = gravity.stopping_conditions.collision_detection
    stellar_evolution, se_light = new_stellar_evolution_codes(stars)
    
    all_se_particles = ParticlesSuperset([se_light.particles, stellar_evolution.particles])
    from_stellar_evolution_to_gravity = all_se_particles.new_channel_to(gravity.particles)
    from_stellar_evolution_to_gravity.copy_attributes(["mass", "radius"])

    print 'Initial conditions:'
    print gravity.particles
    
    mergers = []
    se_merge_products = stellar_evolution.particles[:0]
    
    while gravity.model_time < time_end:
        target_time = min(gravity.model_time + delta_t, time_end)
        gravity.evolve_model(target_time)
        gravity.synchronize_model()
        print "\n   Gravity done, {0} collision(s) detected".format(len(collision_detection.particles(0)))
        try:
            stellar_evolution.evolve_model(gravity.model_time)
        except AmuseException:
            print "This simulations seems to end now... (stellar_evolution returned exception). Mergers so far:"
            print mergers
            raise
        se_light.evolve_model(gravity.model_time)
        print "   Stellar evolution done"
        from_stellar_evolution_to_gravity.copy_attributes(["mass", "radius"])
        
        if collision_detection.is_set():
            gd_primaries = collision_detection.particles(0)
            gd_secondaries = collision_detection.particles(1)
            gd_colliders = gd_primaries + gd_secondaries
            print "   Colliding gravity particles:"
            print gd_colliders
            
            # Are there (low-mass) stars missing in the stellar evolution code (i.e. stored in the 'light' code)?
            missing = gd_colliders.get_intersecting_subset_in(se_light.particles)
            if len(missing) > 0:
                move_from_light_to_real_stellar_evolution(missing, se_light, stellar_evolution)
            
            se_merged, gd_merged = stellar_evolution.merge_colliding(gd_primaries, gd_secondaries, 
                MakeMeAMassiveStar, dict(redirection='file', redirect_file=os.path.join('logs', 'star_collider_code_out.log')), 
                dict(dump_mixed_flag = True, target_n_shells_mixing = 2000),
                return_merge_products=["se", "gd"])
            
            se_merge_products += se_merged
            mergers.extend([[m1, m2, m] for m1, m2, m in zip(gd_primaries.mass.as_quantity_in(units.MSun), gd_secondaries.mass.as_quantity_in(units.MSun), gd_merged.mass)])
            gravity.particles.remove_particles(gd_colliders)
            gravity.particles.add_particles(gd_merged)
            print "   Stellar evolution merge products so far:"
            print se_merge_products
        
        print "   Evolved system (currently {0} particles) up to: {1}".format(len(gravity.particles), gravity.model_time.as_quantity_in(units.Myr))
        max_radius = max(stellar_evolution.particles.radius).as_quantity_in(units.RSun)
        giant = stellar_evolution.particles.select(lambda x : x == max_radius, ["radius"])
        giant_in_gd = giant.get_intersecting_subset_in(gravity.particles)[0]
        print "   Largest star has radius: {0}, and distance to nearest neighbour: {1}".format(max_radius, 
            min(((gravity.particles-giant_in_gd).position - giant_in_gd.position).lengths_squared()).sqrt().as_quantity_in(units.RSun))
        E_kin = gravity.kinetic_energy
        E_pot = gravity.potential_energy
        print "   Total energy K+V: {0}, ratio K/V: {1}".format(E_kin + E_pot, E_kin / E_pot)
        
        write_set_to_file(gravity.particles.copy_to_memory().savepoint(gravity.model_time), 'pop_iii_cluster_N{0}.hdf5'.format(number_of_stars), 'hdf5')

    gravity.stop()
    stellar_evolution.stop()
    se_light.stop()
    
    print "\nSimulation finished. Mergers so far:"
    print mergers
    
if __name__ == '__main__':
    new_working_directory()
    main()
