import os.path
import time
import numpy
from numpy import sin, cos
import cPickle

from amuse.units import nbody_system, units
from amuse.units.quantities import zero
from amuse.datamodel import Particles, ParticlesOverlay
from amuse.io import write_set_to_file, read_set_from_file
from amuse.couple.bridge import Bridge


class GravityHydroStellar(object):
    
    def __init__(self,
            gravity, hydro, stellar,
            gravity_to_hydro, hydro_to_gravity,
            time_step_bridge, 
            time_step_feedback,
            feedback_gas_particle_mass,
            feedback_efficiency = 0.01,
            feedback_radius = 0.01 | units.parsec,
            verbose = True
    ):
        self.gravity = gravity
        self.hydro = hydro
        self.stellar = stellar
        self.gravity_to_hydro = gravity_to_hydro
        self.hydro_to_gravity = hydro_to_gravity
        
        self.time_step_bridge = time_step_bridge
        self.time_step_feedback = time_step_feedback
        self.feedback_gas_particle_mass = feedback_gas_particle_mass
        
        self.feedback_efficiency = feedback_efficiency
        self.feedback_radius = feedback_radius
        self.verbose = verbose
        
        self.bridge = Bridge(verbose=False, timestep=self.time_step_bridge, use_threading=False)
        self.bridge.add_system(self.hydro, (self.gravity_to_hydro,))
        self.bridge.add_system(self.gravity, (self.hydro_to_gravity,))
        
        self.star_particles = self.stars_with_mechanical_luminosity(self.stellar.particles)
        self.total_feedback_energy = zero
        self.current_time = 0.0 | units.Myr
    
    def stars_with_mechanical_luminosity(self, particles):
        result = ParticlesOverlay(particles)
        
        def lmech_function(temperature, mass_loss):
            t4 = numpy.log10(temperature.value_in(units.K)) - 4.0
            v_terminal = (30 + 4000 * t4.clip(0.0, 1.0)) | units.km / units.s
            return 0.5 * mass_loss * v_terminal**2
        
        result.add_calculated_attribute("mechanical_luminosity", 
            lmech_function, attributes_names=["temperature", "mass_loss_wind"])
        result.L_mech = result.mechanical_luminosity
        result.E_mech = (0.0|units.Myr) * result.L_mech
        result.E_mech_last_feedback = result.E_mech
        result.previous_mass = result.mass # Store the mass of the star at the last moment of feedback
        return result
    
    def evolve_model(self, end_time):
        while self.current_time < end_time - 0.5 * self.time_step_feedback:
            self.current_time += self.time_step_feedback
            if self.verbose:
                time_begin = time.time()
                print "GravityHydroStellar: Start evolving..."
            self.bridge.evolve_model(self.current_time)
            self.stellar.evolve_model(self.current_time)
            if self.verbose:
                print "GravityHydroStellar: Evolved to:", self.current_time
                model_times = [getattr(self, code).model_time.as_quantity_in(units.Myr) for code in ["stellar", "bridge", "gravity", "hydro"]]
                print "   (Stellar: {0}, Bridge: {1}, Gravity: {2}, Hydro: {3})".format(*model_times)
                print "GravityHydroStellar: Call mechanical_feedback"
            self.mechanical_feedback(self.time_step_feedback)
            self.store_system_state()
            if self.verbose:
                print "GravityHydroStellar: Iteration took {0} seconds.".format(time.time() - time_begin)
    
    def mechanical_feedback(self, time_step):
        L_mech_new = self.star_particles.mechanical_luminosity
        self.star_particles.E_mech += time_step * 0.5 * (self.star_particles.L_mech + L_mech_new)
        self.star_particles.L_mech = L_mech_new
        self.star_particles.dmass = self.star_particles.previous_mass - self.star_particles.mass
        self.star_particles.n_feedback_particles = numpy.array(self.star_particles.dmass / self.feedback_gas_particle_mass).astype(int)
        
        losers = self.star_particles.select_array(lambda x: x > 0, ["n_feedback_particles"])
        if self.verbose:
            print "GravityHydroStellar: Number of particles providing mechanical feedback during this step:", len(losers)
        if len(losers) == 0:
            return
        channel = self.gravity.particles.new_channel_to(losers)
        channel.copy_attributes(["x","y","z","vx","vy","vz"])
        number_of_new_particles = losers.n_feedback_particles.sum()
        new_sph_all = Particles(number_of_new_particles)
        new_sph_all.mass = self.feedback_gas_particle_mass
        new_sph_all.h_smooth = 0.0 | units.parsec
        offsets = self.draw_random_position_offsets(number_of_new_particles)
        
        i = 0
        for loser in losers:
            i_next = i + loser.n_feedback_particles
            new = new_sph_all[i:i_next]
            new.position = loser.position + offsets[i:i_next]
            new.velocity = loser.velocity
            new.u = self.feedback_efficiency * (loser.E_mech - loser.E_mech_last_feedback) / loser.dmass    
            i = i_next
        
        losers.previous_mass -= self.feedback_gas_particle_mass * losers.n_feedback_particles
        losers.E_mech_last_feedback = losers.E_mech
        channel = losers.new_channel_to(self.gravity.particles)
        channel.copy_attribute("previous_mass", "mass")
        if self.verbose:
            print "GravityHydroStellar: New SPH particles:"
            print new_sph_all
        self.hydro.gas_particles.add_particles(new_sph_all)
        self.total_feedback_energy += (new_sph_all.mass * new_sph_all.u).sum()
    
    def draw_random_position_offsets(self, number_of_new_particles):
        r = numpy.random.uniform(0.0, 1.0, number_of_new_particles)
        theta = numpy.random.uniform(0.0, numpy.pi, number_of_new_particles)
        phi = numpy.random.uniform(0.0, 2*numpy.pi, number_of_new_particles)
        return self.feedback_radius * numpy.array((r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta))).transpose()
    
    @property
    def kinetic_energy(self):
        return self.bridge.kinetic_energy
    
    @property     
    def potential_energy(self):
        return self.bridge.potential_energy
    
    @property
    def thermal_energy(self):
        return self.bridge.thermal_energy
    
    @property
    def feedback_energy(self):
        return self.total_feedback_energy
    
    @property
    def model_time(self):
        return self.bridge.model_time
    
    @property
    def particles(self):
        return self.bridge.particles
    
    @property
    def gas_particles(self):
        return self.bridge.gas_particles
    
    def store_system_state(self):
        snapshot_number = int(0.5 + self.current_time / self.time_step_feedback)
        filename = os.path.join("snapshots", "cluster_snapshot_{0:=06}_".format(snapshot_number))
        stars = Particles(keys=self.star_particles.key)
        self.star_particles.copy_values_of_attributes_to(["mass", "temperature", 
            "stellar_type", "radius", "luminosity", "age", "L_mech", "E_mech", 
            "E_mech_last_feedback", "previous_mass"], stars)
        self.gravity.particles.copy_values_of_attributes_to(["x","y","z","vx","vy","vz"], stars)
        write_set_to_file(stars, filename+"stars.amuse", "amuse")
        write_set_to_file(self.hydro.gas_particles, filename+"gas.amuse", "amuse")
        with open(filename+"info.pkl", "wb") as outfile:
            cPickle.dump(
                (self.gravity.unit_converter, self.hydro.unit_converter,
                self.time_step_bridge, self.time_step_feedback,
                self.feedback_gas_particle_mass,
                self.feedback_efficiency, self.feedback_radius,
                self.verbose, self.total_feedback_energy, self.current_time),
            outfile)
    
    def load_system_state(self, info_file, stars_file):
        stars = read_set_from_file(stars_file, 'amuse')
        channel = stars.new_channel_to(self.star_particles)
        channel.copy_attributes(["L_mech", "E_mech", "E_mech_last_feedback", "previous_mass"])
        
        with open(info_file, "rb") as infile:
            (tmp1, tmp2,
                self.time_step_bridge, self.time_step_feedback,
                self.feedback_gas_particle_mass,
                self.feedback_efficiency, self.feedback_radius,
                self.verbose, self.total_feedback_energy, self.current_time
            ) = cPickle.load(infile)
    
    def stop(self):
        self.bridge.stop()
        self.stellar.stop()
    

