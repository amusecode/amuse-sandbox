import numpy 
import time

from amuse.units import nbody_system, units
from amuse.datamodel import Particle, ParticlesSuperset
from amuse.support.exceptions import AmuseException

from amuse.community.ph4.interface import ph4
from amuse.community.mesa.interface import MESA
from amuse.community.sse.interface import SSE
from amuse.community.mmams.interface import MakeMeAMassiveStar

from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.flatimf import new_flat_mass_distribution


def new_gravity_code(stars, convert_nbody):
    gravity = ph4(convert_nbody, mode='gpu', redirection='file', redirect_file='gravity_code_out.log')
    gravity.stopping_conditions.collision_detection.enable()
    gravity.particles.add_particles(stars)
    return gravity

def new_collision_code():
    star_collider = MakeMeAMassiveStar(redirection='file', redirect_file='star_collider_code_out.log')
    star_collider.initialize_code()
    star_collider.parameters.dump_mixed_flag = True
    star_collider.parameters.target_n_shells_mixing = 2000
    star_collider.commit_parameters()
    return star_collider

def new_stellar_evolution_codes(stars):
    n_low_mass = (stars.mass < 10|units.MSun).sum()
    print n_low_mass, "particles in SSE,", len(stars)-n_low_mass, "particles in MESA"
    
    stellar_evolution = MESA(redirection='file', redirect_file='stellar_evolution_code_out.log', channel_type='ibis', hostname='VU-MC', output_data_root_directory='/home/niels/amuse-data', input_data_root_directory='/home/niels/amuse/data', default_path_to_MESA_data='/home/niels/amuse/src/amuse/community/mesa/src/mesa/data')
    stellar_evolution.parameters.RGB_wind_scheme = 4 # 'Dutch' wind
    stellar_evolution.parameters.AGB_wind_scheme = 4 # 'Dutch' wind
#    stellar_evolution.parameters.dutch_wind_efficiency = 0.8
    stellar_evolution.parameters.metallicity = 0.0
    stellar_evolution.particles.add_particles(stars[n_low_mass:])
    
    se_light = SSE(redirection='file', redirect_file='low_mass_stellar_evolution_code_out.log')
    se_light.parameters.metallicity = 0.0001
    se_light.particles.add_particles(stars[:n_low_mass])
    return stellar_evolution, se_light

def main():
    number_of_stars = 1000
    time_end = 10.0 | units.Myr
    delta_t = 0.01 | units.Myr
    wallclock_begin = time.time()

    numpy.random.seed(12345)
#    masses = new_salpeter_mass_distribution(number_of_stars, 
    masses = new_flat_mass_distribution(number_of_stars, 
        mass_min = 0.1 | units.MSun, mass_max = 100 | units.MSun).sorted()
    
    convert_nbody = nbody_system.nbody_to_si(masses.sum(), 0.01 | units.parsec)
    stars = new_plummer_model(number_of_stars, convert_nbody)
    stars.mass = masses
    stars.scale_to_standard(convert_nbody=convert_nbody)

    gravity = new_gravity_code(stars, convert_nbody)
    star_collider = new_collision_code()
    stellar_evolution, se_light = new_stellar_evolution_codes(stars)
    
    all_se_particles = ParticlesSuperset([se_light.particles, stellar_evolution.particles])
    from_stellar_evolution_to_gravity = all_se_particles.new_channel_to(gravity.particles)
    from_stellar_evolution_to_gravity.copy_attributes(["mass", "radius"])

    print 'ICs:'
    print gravity.particles
    
    mergers = []
    
    while gravity.model_time < time_end:
	wallclock_iteration_begin = time.time()
        target_time = min(gravity.model_time + delta_t, time_end)
	print "Evolving gravity model"
        gravity.evolve_model(target_time)
        gravity.synchronize_model()
        print "\nGravity model evolve done. Collision detected?", gravity.stopping_conditions.collision_detection.is_set()
        try:
	    print "Evolving stellar evolution model"
            stellar_evolution.evolve_model(gravity.model_time)
	    print "Stellar evolution model done"
        except AmuseException:
            print "This simulations seems to end now... (stellar_evolution returned exception). Mergers so far:"
            print mergers
            raise
	print "Evolving low-mass stellar evolution model"
        se_light.evolve_model(gravity.model_time)
        from_stellar_evolution_to_gravity.copy_attributes(["mass", "radius"])
	print "Low-mass stellar evolution model done"
        
        if gravity.stopping_conditions.collision_detection.is_set():
            print '(' + str(len(gravity.stopping_conditions.collision_detection.particles(0))), 'collision(s) detected)'
            gd_primaries = gravity.stopping_conditions.collision_detection.particles(0)
            gd_secondaries = gravity.stopping_conditions.collision_detection.particles(1)
            print gd_primaries + gd_secondaries
            col_primaries = star_collider.particles.add_particles(gd_primaries)
            col_secondaries = star_collider.particles.add_particles(gd_secondaries)
            print "   Distance(s) between colliders:", (gd_primaries.position - gd_secondaries.position).lengths()
            print "   Sum(s) of radii of colliders:", gd_primaries.radius + gd_secondaries.radius
            
            se_colliders = star_collider.particles.get_intersecting_subset_in(stellar_evolution.particles)
            if len(star_collider.particles) != len(se_colliders):
                missing = star_collider.particles.get_intersecting_subset_in(se_light.particles)
                print "   Move the missing particles from the 'light' code to the 'real' stellar evolution code:\n", missing
                vivified = stellar_evolution.particles.add_particles(missing)
                for star in vivified:
                    star.evolve_for(stellar_evolution.model_time)
                print vivified
                se_light.particles.remove_particles(missing)
                se_colliders = star_collider.particles.get_intersecting_subset_in(stellar_evolution.particles)
                if len(star_collider.particles) != len(se_colliders):
                    raise AmuseException("Number of colliding particles does not match the number of particles in", 
                        star_collider.__class__.__name__, "!")
            
            for col_particle, se_particle in zip(star_collider.particles, se_colliders):
                number_of_zones     = se_particle.get_number_of_zones()
                mm1                 = se_particle.get_mass_profile(number_of_zones = number_of_zones)* se_particle.mass
                mass_profile        = se_particle.get_cumulative_mass_profile(number_of_zones = number_of_zones) * se_particle.mass
                density_profile     = se_particle.get_density_profile(number_of_zones = number_of_zones)
                radius_profile      = se_particle.get_radius_profile(number_of_zones = number_of_zones)
                temperature_profile = se_particle.get_temperature_profile(number_of_zones = number_of_zones)
                lum                 = se_particle.get_luminosity_profile(number_of_zones = number_of_zones)
                pressure_profile    = se_particle.get_pressure_profile(number_of_zones = number_of_zones)
                mu_profile          = se_particle.get_mu_profile(number_of_zones = number_of_zones)
                composition_profile = se_particle.get_chemical_abundance_profiles(number_of_zones = number_of_zones)
                
                col_particle.add_shell(mm1,mass_profile, radius_profile, density_profile, 
                    pressure_profile, temperature_profile,lum, mu_profile, composition_profile[0], 
                    composition_profile[1]+composition_profile[2], composition_profile[3], 
                    composition_profile[4], composition_profile[5], composition_profile[6], 
                    composition_profile[7], composition_profile[7]*0.0, composition_profile[7]*0.0)
            
            print "   Copying stellar structure of each particle to", star_collider.__class__.__name__, "done."
            
            for col_primary, col_secondary in zip(col_primaries, col_secondaries):
                merge_product = Particle()
                merge_product.primary = col_primary
                merge_product.secondary = col_secondary
                print "   Merging colliding stars in", star_collider.__class__.__name__
                try:
                    mmams_new_particle = star_collider.merge_products.add_particle(merge_product)
                except AmuseException:
                    print "NOTE: shock_heating_flag = 0"
                    star_collider.set_do_shock_heating_flag(0)
                    mmams_new_particle = star_collider.merge_products.add_particle(merge_product)
                
                stellar_model = mmams_new_particle.internal_structure()
                merge_product.mass = stellar_model.mass[-1] 
                merge_product.radius= stellar_model.radius[-1]
                gd_colliders = (col_primary + col_secondary).get_intersecting_subset_in(gravity.particles)
                merge_product.position = gd_colliders.center_of_mass()
                merge_product.velocity = gd_colliders.center_of_mass_velocity()
                print "   Merge product:\n", merge_product
                mergers.append((gd_colliders.mass, merge_product.mass))
                
                stellar_evolution.new_particle_from_model(stellar_model, 0. | units.Myr, key=merge_product.key)
                gravity.particles.add_particle(merge_product)
                print "   Merge product has been added to gravity and stellar_evolution codes"

            stellar_evolution.particles.remove_particles(star_collider.native_stars)
            gravity.particles.remove_particles(star_collider.native_stars)
            star_collider.particles.remove_particles(star_collider.particles)
            print "   Original colliding particles have been removed from codes"
            
            print "   Stellar evolution particles:"
            print stellar_evolution.particles
            print "   Gravity particles:"
            print gravity.particles

	wallclock_iteration_end = time.time()
        
        print "   Evolved system up to:", gravity.model_time.as_quantity_in(units.Myr)
        print "   System currently contains", len(gravity.particles), "particles."
        max_radius = max(stellar_evolution.particles.radius).as_quantity_in(units.RSun)
        giant = stellar_evolution.particles.select(lambda x : x == max_radius, ["radius"])
        giant_in_gd = giant.get_intersecting_subset_in(gravity.particles)[0]
        print "   Largest star has radius:", max_radius
        print "   Nearest neighbour:", min(((gravity.particles-giant_in_gd).position - giant_in_gd.position).lengths_squared()).sqrt().as_quantity_in(units.RSun)
        E_kin = gravity.particles.kinetic_energy()
        E_pot = gravity.particles.potential_energy()
        print "   Total energy:", E_kin + E_pot
        print "   Energy ratio:", E_kin / E_pot
        print "   Step took", (wallclock_iteration_end - wallclock_iteration_begin) | units.s
        print "   Total time", (wallclock_iteration_end - wallclock_begin) | units.s
	print "   Myr per hour", gravity.model_time.value_in(units.Myr) / ((wallclock_iteration_end - wallclock_begin)/3600)

    print 'Final model:'
    print "   Stellar evolution particles:"
    print stellar_evolution.particles
    print "   Gravity particles:"
    print gravity.particles
    star_collider.stop()
    gravity.stop()
    stellar_evolution.stop()
    se_light.stop()
    
if __name__ == '__main__':
    main()
