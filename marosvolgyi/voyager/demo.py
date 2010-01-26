import sys

from amuse.legacy.bhtree.interface import BHTreeInterface, BHTree
from amuse.legacy.hermite0.interface import HermiteInterface, Hermite
from amuse.legacy.phiGRAPE.interface import PhiGRAPEInterface, PhiGRAPE
import amuse.legacy.twobody.twobody as twobody

from sandbox.marosvolgyi.viewer import planetarium
from amuse.support.io.horizons import LoadStar, NewStar

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units

from pylab import *

from datetime import date, timedelta

class Horizons():
    def __init__(self):

        self.Sun = LoadStar('sun')
        self.Sun.mass = units.MSun(1.0)
        self.Sun.radius = units.RSun(1.0)

        self.Mercury = LoadStar('mercury')
        self.Mercury.mass = units.kg( 3.302e23)
        self.Mercury.radius = units.km(2440) 

        self.Venus = LoadStar('venus')
        self.Venus.mass = units.kg(48.685e23)
        self.Venus.radius = units.km(6051.8) 
        
        self.Earth = LoadStar('earth')
        self.Earth.mass = units.kg(5.9736e24)
        self.Earth.radius = units.km(6371) 

        self.Moon = LoadStar('moon')
        self.Moon.mass = units.kg(7.35e22)
        self.Moon.radius = units.km(1738) 

        self.Mars = LoadStar('mars')
        self.Mars.mass = units.kg(6.4185e23)
        self.Mars.radius = units.km(3389.9)

        self.Jupiter = LoadStar('jupiter')
        self.Jupiter.mass = units.kg(1898.13e24)
        self.Jupiter.radius = units.km(71492)
        #self.Jupiter.radius = units.km(0.000001)

        self.Saturn = LoadStar('saturn')
        self.Saturn.mass = units.kg(5.68319e26 )
        self.Saturn.radius = units.km(58232)

        self.Uranus = LoadStar('uranus')
        self.Uranus.mass = units.kg(86.8103e24 )
        self.Uranus.radius = units.km(26000)

        self.Neptune = LoadStar('neptune')
        self.Neptune.mass = units.kg(102.41e24 )
        self.Neptune.radius = units.km(25000)

        self.VoyagerI = LoadStar('voyagerI')
        self.VoyagerI.mass = units.kg(500)
        self.VoyagerI.radius = units.m(20)

        self.VoyagerII = LoadStar('voyagerII')
        self.VoyagerII.mass = units.kg(1000)
        self.VoyagerII.radius = units.m(20)
    

class SolarSystemModel(object):

    def __init__(self):
        self.load_integrator_and_units()
        self.planets = []
        self.planets_data = Horizons()

        #Notable events
        self.start_date = date(1971, 10, 26)
        self.voyagerI_launch_date = date(1977, 9, 7)
        self.stop_date = date(1989, 10, 26)

        self.model_t0 = self.start_date

        self.planets = core.Stars(10)
        self.voyagerI = core.Particle()

        #set I.C. using Horizons database
        self.set_IC_at_date([self.planets_data.Sun, self.planets_data.Mercury, self.planets_data.Venus, 
                             self.planets_data.Earth, self.planets_data.Moon, self.planets_data.Mars, 
                             self.planets_data.Jupiter, self.planets_data.Saturn, self.planets_data.Uranus,
                             self.planets_data.Neptune], 
                            self.planets,
                            self.start_date)

        self.set_IC_at_date([self.planets_data.VoyagerI],self.voyagerI.as_set(), self.voyagerI_launch_date)

        self.planets.synchronize_to(self.code.particles)

        self.planets_from_code_to_model = self.code.particles.new_channel_to(self.planets)

        #later we will add particles
        self.all_particles = self.planets.copy()

        self.P = planetarium.SolarSystemView((1600,1000))

    def load_integrator_and_units(self):
        self.convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        #self.code = PhiGRAPE(self.convert_nbody)
        #self.code.parameters.epsilon_squared = 0.0 | units.AU**2
        #self.code.set_eta(0.01,0.02)
        #self.code.dt_dia = 5000

        #self.code = BHTree(self.convert_nbody)
        #self.code.parameters.epsilon_squared = 0.000000001 | units.AU**2
        #self.instance.set_eta(0.0001,0.0002)
        #self.code.dt_dia = 10
        #self.code.timestep = 0.002 
        #instance.use_self_gravity = 0

        self.code = Hermite(self.convert_nbody)
        self.code.parameters.epsilon_squared = 0.0 | units.AU**2
        self.code.dt_dia = 5000

        self.code.setup_module()

    def set_IC_at_date(self, data, particles, IC_date):

        for i, body in enumerate(data):

            particles[i].mass = body.mass
            particles[i].radius = body.radius
            r, v = body.get_vectors_at_date(IC_date)
            particles[i].position = units.AU([r[0],r[1],r[2]])
            particles[i].velocity = units.AUd([v[0],v[1],v[2]])

    def approach(self, voyagerI, chkdistto):
        dR = voyagerI.position - chkdistto.position

        if (dot(dR,dR)**0.5) < 0.0120|units.AU:
            return True
        else:
            return False

    def retreat(self, state):
        dist = (state['x']**2+state['y']**2+state['z']**2)**0.5
        if (dist|units.m).value_in(units.AU) > 0.0120:
            return True
        else:
            return False


    def evolve(self, start_date, end_date, step):

        for x in arange(start_date, end_date, step):
            self.code.evolve_model((x-self.model_t0.toordinal())|units.day)
            if not self.code.get_indices_of_colliding_particles()['index_of_particle1'] == -1:
                pass
                
            #self.code.update_particles(self.planets)
            #self.planets.savepoint(x)

    def other_evolve(self, t):
        pass

    def remove_pair(self, a_planet):

        self.code.particles.remove_particle(self.voyagerI)
        self.code.particles.remove_particle(a_planet) 
        self.code.particles.synchronize_to(self.all_particles)
        self.code.update_particles(self.all_particles)
        #planets are heavy wtr spacecraft 
        #so we put the planet back instead of the c.o.g. sys
        self.all_particles.add_particle(a_planet)
        self.all_particles.synchronize_to(self.code.particles)
        self.a_planet_from_code_to_model =  self.code.particles.new_channel_to(a_planet.as_set())
        self.code.update_particles(self.all_particles)

    def restore_pair(self, state, a_planet):
        self.voyagerI = core.Particle()
        voyagerI_attrib = self.voyagerI.as_set()
        rel_vr = ([state['vx'],state['vy'],state['vz']]|units.m/units.s).as_quantity_in(units.AUd)
        rel_r = ([state['x'],state['y'],state['z']]|units.m).as_quantity_in(units.AU)
        voyagerI_attrib[0].velocity = a_planet.velocity + rel_vr
        voyagerI_attrib[0].position = a_planet.position + rel_r
        
        self.all_particles.add_particle(self.voyagerI)
        self.all_particles.synchronize_to(self.code.particles)
        self.voyagerI_from_code_to_model = self.code.particles.new_channel_to(self.voyagerI.as_set())
    
    def setup_two_body_system(self, a_planet):
        m = 1898.13e24#kg
        #ra = 71492000.0#m
        ra = 1.0#m

        dr = self.voyagerI.position - a_planet.position
        dv = self.voyagerI.velocity - a_planet.velocity
        del self.voyagerI
        tb = twobody.twobody()
        tb.new_particle(m, 
                        ra, 
                        dr[0].value_in(units.m),
                        dr[1].value_in(units.m),
                        dr[2].value_in(units.m),
                        dv[0].value_in(units.m/units.s)*1.5,
                        dv[1].value_in(units.m/units.s)*1.5,
                        dv[2].value_in(units.m/units.s)*1.5)
        return tb

    def handle_viewer(self):
        self.P.renderamuse(self.all_particles)
        self.P.handle_events()
        if not self.P.go:
            self.P.go = True
            return False
        else:
            return True

    def add_voyagers_with_channel(self):
        self.all_particles.add_particle(self.voyagerI)
        self.all_particles.synchronize_to(self.code.particles)
        self.voyagerI_from_code_to_model = I.code.particles.new_channel_to(self.voyagerI.as_set())

    def run_stage_1(self):
        
        for t in range(self.start_date.toordinal(), self.voyagerI_launch_date.toordinal(), 1):

            self.evolve(t,t+1,1)
            self.code.update_particles(self.all_particles)        
            
            self.planets_from_code_to_model.copy()
            
            if not self.handle_viewer():
                break

        s = raw_input()

    def run_stage_2(self):

        self.add_voyagers_with_channel()

        two_body_mode = False
        timestep = .2

        for t in arange(self.voyagerI_launch_date.toordinal(), self.stop_date.toordinal(), timestep):

            self.evolve(t,t+timestep,timestep)
            self.code.update_particles(self.all_particles)        

            self.planets_from_code_to_model.copy()

            
            if not two_body_mode:

                self.voyagerI_from_code_to_model.copy()

                for a_planet in self.planets:
                    
                    if self.approach(self.voyagerI, a_planet):
                        print "approach"
                        self.remove_pair(a_planet)
                        two_body_mode = True
                        tb = self.setup_two_body_system(a_planet)
                        days = 1
                        break

            elif two_body_mode:

                tb.evolve(3600.0*24*1*days)
                days += timestep
                self.a_planet_from_code_to_model.copy()
                state, err = tb.get_state(0)

                if self.retreat(state):
                    two_body_mode = False
                    print "retreat"
                    self.restore_pair(state, a_planet)
                    timestep = .2
                    del tb
                    
            if not self.handle_viewer():
                break

        s = raw_input()

        
    def cleanup(self):
        del self.instance



if __name__ == '__main__':

    I = SolarSystemModel()
    I.run_stage_1()
    I.run_stage_2()
