import numpy 

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
    gravity.initialize_code()
    gravity.parameters.set_defaults()
    gravity.parameters.epsilon_squared = 0.0 | nbody_system.length**2
    gravity.commit_parameters()
    gravity.stopping_conditions.collision_detection.enable()
    gravity.particles.add_particles(stars)
    gravity.commit_particles()
    return gravity

def main():
    number_of_stars = 900
    numpy.random.seed(12345)
#    masses = new_salpeter_mass_distribution(number_of_stars, 
    masses = new_flat_mass_distribution(number_of_stars, 
        mass_min = 0.1 | units.MSun, mass_max = 100 | units.MSun).sorted()
    
    convert_nbody = nbody_system.nbody_to_si(masses.sum(), 0.01 | units.parsec)
    stars = new_plummer_model(number_of_stars, convert_nbody)
    stars.mass = masses
    stars.scale_to_standard(convert_nbody=convert_nbody)

    gravity = new_gravity_code(stars, convert_nbody)
        
    star_collider = MakeMeAMassiveStar(redirection='file', redirect_file='star_collider_code_out.log')
    star_collider.initialize_code()
    star_collider.parameters.dump_mixed_flag = True
    star_collider.parameters.target_n_shells_mixing = 2000
    star_collider.commit_parameters()
    
    n_low_mass = len(numpy.where(masses < 10|units.MSun)[0])
    se_light = SSE(redirection='file', redirect_file='low_mass_stellar_evolution_code_out.log')
    se_light.parameters.metallicity = 0.0001
    se_light.particles.add_particles(stars[:n_low_mass])
    
    stellar_evolution = MESA(redirection='file', redirect_file='stellar_evolution_code_out.log')
    stellar_evolution.parameters.RGB_wind_scheme = 4 # 'Dutch' wind
    stellar_evolution.parameters.AGB_wind_scheme = 4 # 'Dutch' wind
#    stellar_evolution.parameters.dutch_wind_efficiency = 0.8
    stellar_evolution.parameters.metallicity = 0.0
    stellar_evolution.particles.add_particles(stars[n_low_mass:])
    
    all_se_particles = ParticlesSuperset([se_light.particles, stellar_evolution.particles])
    
    from_stellar_evolution_to_gravity = stellar_evolution.particles.new_channel_to(gravity.particles)
    from_stellar_evolution_to_gravity.copy_attributes(["mass", "radius"])
    
    from_light_se_to_gravity = se_light.particles.new_channel_to(gravity.particles)
    from_light_se_to_gravity.copy_attributes(["mass", "radius"])

    time_end = 10.0 | units.Myr
    delta_t = 0.01 | units.Myr

    print 'ICs:'
    print gravity.particles
    
    mergers = []
    
    while gravity.model_time < time_end:
        target_time = min(gravity.model_time + delta_t, time_end)
        gravity.evolve_model(target_time)
        gravity.synchronize_model()
        print "\nCollision detected?", gravity.stopping_conditions.collision_detection.is_set(),
        try:
            stellar_evolution.evolve_model(gravity.model_time)
        except AmuseException:
            print "This simulations seems to end now... (stellar_evolution returned exception). Mergers so far:"
            print mergers
            raise
        se_light.evolve_model(gravity.model_time)
        from_stellar_evolution_to_gravity.copy_attributes(["mass", "radius"])
        from_light_se_to_gravity.copy_attributes(["mass", "radius"])
        
        # JUST FOR TESTING!!!
#        gravity.particles.radius *= 1000
        
        if not gravity.stopping_conditions.collision_detection.is_set():
            print '(no collision detected)'
        else:
            print '(' + str(len(gravity.stopping_conditions.collision_detection.particles(0))), 'collision(s) detected)'
            time_col = gravity.model_time
            
            print gravity.particles
            gd_primaries = gravity.stopping_conditions.collision_detection.particles(0)
            gd_secondaries = gravity.stopping_conditions.collision_detection.particles(1)
            print gd_primaries.key, gd_secondaries.key
            col_primaries = star_collider.particles.add_particles(gd_primaries)
            col_secondaries = star_collider.particles.add_particles(gd_secondaries)
            print "   Distance(s) between colliders:", (gd_primaries.position - gd_secondaries.position).lengths()
            print "   Sum(s) of radii of colliders:", gd_primaries.radius + gd_secondaries.radius
            
            # Store info from non-interacting particles to check they're still there after merger
            if len(all_se_particles) > 2:
                non_colliding_se_masses = (all_se_particles - col_primaries - col_secondaries).mass
                non_colliding_gd_posxs = (gravity.particles - col_primaries - col_secondaries).x
            
            print "   Number of particles in", star_collider.__class__.__name__ + ":", len(star_collider.particles)
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
                species_names       = se_particle.get_names_of_species()
                
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
                
                print "   Merging colliding stars done"
                print "   Number of zones in merge product:", mmams_new_particle.number_of_zones
                
                stellar_model = mmams_new_particle.internal_structure()
#~                print stellar_model.mass[[0, -1]], stellar_model.radius[[0, -1]], stellar_model.temperature[[0, -1]], stellar_model.X_H[[0, -1]]
#~                print "stellar_model:"
                
                merge_product.mass = stellar_model.mass[-1] 
                merge_product.radius= stellar_model.radius[-1]
    #            merge_product.radius= 10000.*stellar_model.radius[-1]
                gd_colliders = (col_primary + col_secondary).get_intersecting_subset_in(gravity.particles)
                print gd_colliders
                merge_product.position = gd_colliders.center_of_mass()
                merge_product.velocity = gd_colliders.center_of_mass_velocity()
                print "   Merge product:\n", merge_product.as_set()
                mergers.append((gd_colliders.mass, merge_product.mass))
                
                try:
                    stellar_evolution.new_particle_from_model(stellar_model, 0. | units.Myr, key=merge_product.key)
                except AmuseException:
                    print 'NOTE: Failing to import model:'
                    stellar_evolution.particles.add_particle(merge_product)
                
                gravity.particles.add_particle(merge_product)
                print "   Merge product has been added to gravity and stellar_evolution codes"

#~            print
#~            print 'After merger the code evole stars:'
#~            print stellar_evolution.particles
#~            print star_collider.particles

#CHECK che sia corretto che particelle vengano rimosse dopo aver aggiunto il merger remnant! Come si comportano gli indici????
            stellar_evolution.particles.remove_particles(star_collider.native_stars)
            gravity.particles.remove_particles(star_collider.native_stars)
            star_collider.particles.remove_particles(star_collider.native_stars)
            if len(gravity.particles) < 2:
                gravity.stopping_conditions.collision_detection.disable()
            print "   Original colliding particles have been removed from codes"
            
            # Check whether the non-interacting particles are still there after merger
            if len(all_se_particles) > 1:
                if abs(non_colliding_se_masses - 
                    (all_se_particles - star_collider.merge_products).mass).sum() > 0 | units.MSun:
                    raise AmuseException("Wrong particles seem to have been removed!")
                if abs(non_colliding_gd_posxs - 
                    (gravity.particles - star_collider.merge_products).x).sum() > 0 | units.m:
                    raise AmuseException("Wrong particles seem to have been removed!")
            star_collider.particles.remove_particles(star_collider.merge_products)
            
            print "   Stellar evolution particles:"
            print stellar_evolution.particles
            print "   Gravity particles:"
            print gravity.particles
        
#~            print gravity.stopping_conditions.collision_detection.is_set()
#Adesso devo rimettere a zero la stopping condition immagino........

        print "   Evolved system up to:", gravity.model_time.as_quantity_in(units.Myr)
        print "   System currently contains", len(gravity.particles), "particles."
        max_radius = max(stellar_evolution.particles.radius).as_quantity_in(units.RSun)
        giant = stellar_evolution.particles.select(lambda x : x == max_radius, ["radius"])
        giant_in_gd = giant.get_intersecting_subset_in(gravity.particles)[0]
        print "   Largest star has radius:", max_radius, giant.radius
        print "   Nearest neighbour:", min(((gravity.particles-giant_in_gd).position - giant_in_gd.position).lengths_squared()).sqrt().as_quantity_in(units.RSun)
        E_kin = gravity.particles.kinetic_energy()
        E_pot = gravity.particles.potential_energy()
        print "   Total energy:", E_kin + E_pot
        print "   Energy ratio:", E_kin / E_pot

#--------------QUI TERMINA CICLO SUI J (LABEL DELLE COLLISIONS!)-------------------------------


    print 'Final model:'
    print "   Stellar evolution particles:"
    print stellar_evolution.particles
    print "   Gravity particles:"
    print gravity.particles
    star_collider.stop()
    gravity.stop()
    stellar_evolution.stop()
    
    
if __name__ == '__main__':
    main()
