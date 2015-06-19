from amuse import datamodel
from amuse.units import nbody_system

from amuse import datamodel

from amuse.units import nbody_system


from collections import namedtuple



import numpy




class DiskCode(object):
    
    def __init__(self):
        self.central_particles = datamodel.Particles()
        self.orbiters = datamodel.Particles()
        self.partciles = datamodel.ParticlesSuperset([self.central_particles, self.orbiters])
        self.begin_time = 0 | nbody_system.time
        self.time = self.begin_time
        self.dtparameter = 0.03
        self.G = nbody_system.G
    

    def evolve_model(self, end_time):
        integrator = HermiteIntegrator(self.central_particles, self.orbiters, self.G)
        integrator.setup(False)
        while self.time < end_time:
            integrator.step(self.dt)
            self.time += self.dt
            self.dt = self.calculcate_timestep_from_collision_time(self.orbiters.collection_attributes.collision_time)
        integrator.update_set()




    def commit_particles(self):
        integrator = HermiteIntegrator(self.central_particles, self.orbiters, self.G)
        integrator.setup(True)
        integrator.calculate_derivatives()
        integrator.update_set()
        self.time = self.begin_time
        self.dt = self.calculcate_timestep_from_collision_time(self.orbiters.collection_attributes.collision_time)
        




    def get_acceleration_jerk_and_collision_time(self, iparticles, jparticles, G = nbody_system.G):
        ipos = iparticles.position
        ivel = iparticles.velocity
        imass = iparticles.mass
        jpos = jparticles.position
        jvel = jparticles.velocity
        jmass = jparticles.mass
        
        return self.get_acc_jerki(imass, ipos, ivel, jmass, jpos, jvel, G)

    def get_time_estimate(self, imass, jmass, dr_squared, dv_squared, acc, G):
        estimate1 = ((dr_squared ** 2) / (dv_squared ** 2)).sum(1)
        total_mass = imass + jmass
        estimate2 = ( dr_squared.sum(1)) / ( G * G * (acc.lengths_squared() * (total_mass * total_mass).sum(1)))
        return min(estimate1.min(), estimate2.min()).sqrt().sqrt()

    def get_acc_jerki(self, imass, ipos, ivel, jmass, jpos, jvel, G = nbody_system.G):
        min_G = -G
        n = len(jpos)
        ni = len(ipos)
        newshape =(n, 1, 3)
        jpos = jpos.reshape(newshape)         
        jvel = jvel.reshape(newshape)         
        jmass = jmass.reshape([n,1])
        dpos = jpos - ipos
        dvel = jvel - ivel
        dr_squared = dpos.lengths_squared()
        dv_squared = dvel.lengths_squared()
        dr = dr_squared**0.5
        dr_qubed = dr * dr_squared
        drdv = (dpos * dvel).sum(2)
        drdv_dr_squared = drdv / dr_squared
        m_div_dr_qubed = (imass / dr_qubed).reshape([n,ni,1])
        acc = min_G * (dpos * m_div_dr_qubed).sum(1)
        jerk = min_G * ( (dvel - 3 * drdv_dr_squared.reshape([n,ni,1]) * dpos) * m_div_dr_qubed).sum(1) ;
        
        
        return acc, jerk,self.get_time_estimate(imass, jmass, dr_squared, dv_squared, (dpos / dr_qubed.reshape([n,ni,1])).sum(1), G)


    def calculcate_timestep_from_collision_time(self, collision_time):
        return self.dtparameter * collision_time

    def calculate_derivatives(self, particles):
        particles.acceleration, particles.jerk, particles.collection_attributes.collision_time = self.get_acceleration_jerk_and_collision_time(self.central_particles, particles, G=self.G)
        


    def predict(self, particles, dt):
        def delta_position(p):
            return (p.velocity * dt) + (p.acceleration*dt*dt/2.0) + (p.jerk*dt*dt*dt/6.0)
        def delta_velocity(p):
            return (p.acceleration * dt) + (p.jerk*dt*dt/2.0);
        
        predicted_set = particles.copy()
        predicted_set.position += delta_position(particles)
        predicted_set.velocity += delta_velocity(particles)
        self.calculate_derivatives(predicted_set)
        return predicted_set

    def correct(self, particles, predicted_particles, dt):
        old_vel = particles.velocity
        particles.velocity += (particles.acceleration + predicted_particles.acceleration) * dt / 2.0
        particles.velocity += (particles.jerk         - predicted_particles.jerk) * dt * dt / 12.0
        
        particles.position += (old_vel                + particles.velocity) * dt / 2.0
        particles.position += (particles.acceleration - predicted_particles.acceleration) * dt * dt / 12.0
        
        particles.acceleration = predicted_particles.acceleration
        particles.jerk         = predicted_particles.jerk
        particles.collection_attributes.collision_time = predicted_particles.collection_attributes.collision_time
    def evolve_model_old(self, end_time):
        while self.time < end_time:
            self.correct(self.orbiters, self.predict(self.orbiters, self.dt), self.dt)
            self.time += self.dt
            self.dt = self.calculcate_timestep_from_collision_time(self.orbiters.collection_attributes.collision_time)




class HermiteIntegrator(object):
    
    def __init__(self, iparticles, jparticles , G):
        self.iparticles = iparticles
        self.jparticles = jparticles
        self.G = G
    

    def setup(self, initial = True):
        self.ipos = self.iparticles.position
        self.ivel = self.iparticles.velocity
        self.imass = self.iparticles.mass
        self.jpos = self.jparticles.position
        self.jvel = self.jparticles.velocity
        n = len(self.jparticles)
        self.jmass = self.jparticles.mass.reshape([n,1])
        
        self.ikey = self.iparticles.key
        self.jkey = self.jparticles.key.reshape([n,1])
        self.mask = self.ikey == self.jkey     # (~(self.ikey == self.jkey)) * 1.0
        
        
        
        self.length_unit = self.jpos.unit
        self.time_unit = self.length_unit / self.jvel.unit                                                                              
        self.mass_unit = self.imass.unit
        self.speed_unit = self.length_unit / self.time_unit
        self.wu = 1
        if self.wu:
            self.ipos = self.ipos.value_in(self.length_unit)
            self.jpos = self.jpos.value_in(self.length_unit)
            self.ivel = self.ivel.value_in(self.speed_unit)
            self.jvel = self.jvel.value_in(self.speed_unit)
            self.imass = self.imass.value_in(self.mass_unit)
            self.jmass = self.jmass.value_in(self.mass_unit)
        if self.wu:
          
            if not initial:
                self.acceleration = self.jparticles.acceleration.value_in(self.length_unit * self.time_unit ** -2)
                self.jerk = self.jparticles.jerk.value_in(self.length_unit * self.time_unit ** -3)
                self.collision_time = self.jparticles.collection_attributes.collision_time.value_in(self.time_unit)
        else:
            if not initial:
                self.acceleration = self.jparticles.acceleration
                self.jerk = self.jparticles.jerk
                self.collision_time = self.jparticles.collection_attributes.collision_time
            
        total_mass = self.imass + self.jmass
        self.total_mass_squared = (total_mass * total_mass).sum(1)
        if self.wu:
            self.min_G = (-self.G).value_in(self.length_unit**3 * self.mass_unit**-1 * self.time_unit**-2)
        else:
            self.min_G = -self.G
        self.G_squared = self.min_G*self.min_G



    def step(self, dt):
        if self.wu:
            dt = dt.value_in(self.time_unit)
        self.store_start_values()
        self.predict(dt)
        self.calculate_derivatives()
        self.correct(dt)




    def get_acceleration_jerk_and_collision_time(self):
        return self.get_acc_jerki(self.imass, self.ipos, self.ivel, self.jmass, self.jpos, self.jvel, self.G)


    def get_time_estimate(self, dr_squared, dv_squared, acc, G):
        tmp = ((dr_squared ** 2) / (dv_squared ** 2))
        tmp[self.mask] = 0
        estimate1 = tmp.sum(1)
        estimate2 = ( dr_squared.sum(1)) / ( self.G_squared * (self.lengths_squared(acc) * self.total_mass_squared  ))
        return numpy.sqrt(numpy.sqrt( min(estimate1.min(), estimate2.min()) ))



    def get_acc_jerki(self, imass, ipos, ivel, jmass, jpos, jvel, G = nbody_system.G):
        n = len(jpos)
        ni = len(ipos)
        newshape =(n, 1, 3)
        mask = self.mask
        jpos = jpos.reshape(newshape)                                                                                                                     
        jvel = jvel.reshape(newshape)                                                                                                         
        dpos = jpos - ipos
        dvel = jvel - ivel
        dr_squared = self.lengths_squared(dpos)
        dv_squared = self.lengths_squared(dvel)
        dr = dr_squared**0.5
        dr_qubed = dr * dr_squared
        drdv = (dpos * dvel).sum(2)
        one_div_dr = 1.0 / dr_qubed
        one_div_dr[mask] = 0
        
        drdv_dr_squared =(drdv / dr_squared)
        drdv_dr_squared[mask] = 0
        
        newshape =(n,ni,1)
        m_div_dr_qubed = imass * one_div_dr
        m_div_dr_qubed = (m_div_dr_qubed).reshape(newshape)
        acc = self.min_G * (dpos * m_div_dr_qubed).sum(1)
        jerk = self.min_G * ( (dvel - 3.0 * drdv_dr_squared.reshape(newshape) * dpos) * m_div_dr_qubed).sum(1) ;
        
        return acc, jerk, self.get_time_estimate(dr_squared, dv_squared, (dpos * (one_div_dr).reshape(newshape)).sum(1), G)




    def calculate_derivatives(self):
        self.acceleration, self.jerk, self.collision_time = self.get_acceleration_jerk_and_collision_time()
        



    def predict(self, dt):
        dt_squared = (dt**2/2.0)
        dt_cubed = (dt**3/6.0)

        self.jpos = self.jpos + (self.jvel * dt) + (self.acceleration*dt_squared) + (self.jerk*dt_cubed)
        self.jvel = self.jvel + (self.acceleration * dt) + (self.jerk*dt_squared)



    def correct(self, dt):
        self.jvel = self.start_jvel + ((self.start_acceleration + self.acceleration) * (dt / 2.0)) + ((self.start_jerk         - self.jerk) * (dt ** 2 / 12.0))
        self.jpos = self.start_jpos + ((self.start_jvel         + self.jvel) * (dt / 2.0))         + ((self.start_acceleration - self.acceleration) * (dt**2 / 12.0))


    def store_start_values(self):
        
        self.start_jpos = self.jpos
        self.start_jvel = self.jvel
        self.start_acceleration = self.acceleration
        self.start_jerk = self.jerk
        


    def update_set(self):
        if self.wu:
            self.jparticles.velocity = self.jvel | self.speed_unit
            self.jparticles.position = self.jpos | self.length_unit
            self.jparticles.acceleration = self.acceleration | (self.length_unit * (self.time_unit ** -2))
            self.jparticles.jerk = self.jerk | (self.length_unit * (self.time_unit ** -3))
            self.jparticles.collection_attributes.collision_time = self.collision_time | self.time_unit
        else:
            self.jparticles.velocity = self.jvel             
            self.jparticles.position = self.jpos
            self.jparticles.acceleration = self.acceleration             
            self.jparticles.jerk = self.jerk
            self.jparticles.collection_attributes.collision_time = self.collision_time



    def lengths_squared(self, value):
        if self.wu:
            return (value*value).sum(value.ndim - 1)
        else:
            return value.lengths_squared()

    def movejtoi(self): 
        self.ipos = self.jpos.copy()
        self.ivel = self.jvel.copy()



class HermiteCode(object):
    
    def __init__(self):
        self.particles = datamodel.Particles()
        self.orbiters = self.particles
        self.central_particles = self.particles
        self.begin_time = 0 | nbody_system.time
        self.time = self.begin_time
        self.dtparameter = 0.03
        self.G = nbody_system.G
    


    def evolve_model(self, end_time):
        integrator = HermiteIntegrator(self.particles, self.particles, self.G)
        integrator.setup(False)
        while self.time < end_time:
            integrator.step(self.dt)
            self.time += self.dt
            self.dt = self.calculcate_timestep_from_collision_time(self.orbiters.collection_attributes.collision_time)
            integrator.movejtoi()
        integrator.update_set()





    def commit_particles(self):
        integrator = HermiteIntegrator(self.central_particles, self.orbiters, self.G)
        integrator.setup(True)
        integrator.calculate_derivatives()
        integrator.update_set()
        self.time = self.begin_time
        self.dt = self.calculcate_timestep_from_collision_time(self.orbiters.collection_attributes.collision_time)
        




    def calculcate_timestep_from_collision_time(self, collision_time):
        return self.dtparameter * collision_time


