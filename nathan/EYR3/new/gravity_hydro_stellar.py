import numpy
import time
from numpy import sin, cos

from amuse.units import nbody_system, units
    
from amuse.datamodel import Particles, ParticlesOverlay
from amuse.units.quantities import zero

from amuse.couple.bridge import Bridge

from amuse.ext.evrard_test import uniform_unit_sphere



class GravityHydroStellar(object):
    
    def __init__(self,
            gravity, hydro, stellar,
            gravity_to_hydro, hydro_to_gravity,
            time_step_bridge, 
            time_step_feedback,
            feedback_gas_particle_mass,
            feedback_efficiency = 0.01,
            feedback_radius = 0.01 | units.parsec,
#~            feedback_safety = 1.0e-4 | 1.0e51 * units.erg,
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
        
        self.bridge = Bridge(verbose=verbose, timestep=self.time_step_bridge, use_threading=False)
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
        
        result.add_calculated_attribute("mechanical_luminosity", lmech_function, attributes_names=["temperature", "mass_loss_wind"])
        result.L_mech = result.mechanical_luminosity
        result.E_mech = (0.0|units.Myr) * result.L_mech
        result.E_mech_last_feedback = result.E_mech
        result.previous_mass = result.mass # Store the mass of the star at the last moment of feedback
        return result
    
    def evolve_model(self, end_time):
        while self.current_time < end_time - 0.5 * self.time_step_feedback:
            self.current_time += self.time_step_feedback
            if self.verbose:
		beginning = time.time()
                print "GravityHydroStellar: Start evolving..."
            self.bridge.evolve_model(self.current_time)
#~            print self.bridge.model_time,',',self.sph.model_time,self.grav.model_time
            self.stellar.evolve_model(self.current_time)
            if self.verbose:
                print "GravityHydroStellar: Evolved to:", self.current_time
                print "GravityHydroStellar: Call mechanical_feedback"
            self.mechanical_feedback(self.time_step_feedback)
	    if self.verbose:
		end = time.time()
		print 'evolve to', self.current_time, 'took:', (end - beginning), 'seconds'
    
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
    
    def dump_system_state(self,filename):    
        from amuse.io import write_set_to_file
        import cPickle
        write_set_to_file(self.grav.particles,filename+".grav","amuse",append_to_file=False)  
        write_set_to_file(self.gas_particles,filename+".gas","amuse",append_to_file=False)  
        write_set_to_file(self.star_particles,filename+".evo","amuse",append_to_file=False)  
        write_set_to_file(self.star_particles_addition,filename+".add","amuse",append_to_file=False)  
        f=open(filename+".data",'wb')
        print self.total_feedback_energy
        cPickle.dump((self.codes,
                      self.conv,
                      self.parameters,
                      self.mgas,
                      self.feedback_efficiency,
                      self.feedback_radius,
                      self.time,
                      self.dt_feedback,
                      self.dt_fast,
                      self.total_feedback_energy.in_(1.e51*units.erg),
                      self.feedback_safety,
                      self.feedback_dt,
                      self.feedback_period,
                      self.feedback_lasttime
                      ),f)
        f.close()

    @classmethod
    def load_system_state(cls,filename,new_gas_options=(),
               grav_code_extra=dict(mode='gpu', redirection='none'),
               gas_code_extra=dict(number_of_workers=3,use_gl=False, redirection='none'),
               se_code_extra=dict(redirection='none'),
               grav_couple_code_extra=dict()):    
        from amuse.io import read_set_from_file
        import cPickle
        star_parts=read_set_from_file(filename+".grav",'amuse')
        gas_parts=read_set_from_file(filename+".gas",'amuse')
        evo=read_set_from_file(filename+".evo",'amuse')
        add=read_set_from_file(filename+".add",'amuse')
        f=open(filename+".data",'r')
        data=cPickle.load(f)
        f.close()
        gas_code=data[0][0]
        grav_code=data[0][1]
        se_code=data[0][2]
        grav_couple_code=data[0][3]
        conv=data[1]
        gravp=data[2][0]
        gasp=data[2][1]
        cp=data[2][2]
        mgas=data[3]
        fe=data[4]
        fr=data[5]
        to=data[6]
        dt_feedback=data[7]
        dt_fast=data[8]
        tfe=data[9]
        fs=data[10]
        fdt=data[11]
        fp=data[12]
        flt=data[13]
        
        gasp=gasp+new_gas_options    
            
        return conv,cls(grav_code,gas_code,se_code,grav_couple_code, 
                   conv,mgas,star_parts,gas_parts,dt_feedback, dt_fast,
                   grav_parameters=gravp,gas_parameters=gasp,couple_parameters=cp,
                   feedback_efficiency=fe, feedback_radius=fr,
                   total_feedback_energy=tfe,evo_particles=evo, star_particles_addition=add,
                   start_time_offset=to,feedback_safety=fs,
                   feedback_dt=fdt, feedback_period=fp,
                   feedback_lasttime=flt,
               grav_code_extra=grav_code_extra,
                   gas_code_extra=gas_code_extra,
                   se_code_extra=se_code_extra,
                   grav_couple_code_extra=grav_couple_code_extra)
    
