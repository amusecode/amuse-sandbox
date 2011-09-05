import sys

from amuse.legacy.bhtree.interface import BHTreeInterface, BHTree
from amuse.legacy.hermite0.interface import HermiteInterface, Hermite
from amuse.legacy.phiGRAPE.interface import PhiGRAPEInterface, PhiGRAPE
import amuse.legacy.twobody.twobody as twobody

from sandbox.marosvolgyi.viewer import planetarium
from amuse.io.horizons import LoadStar, NewStar

from amuse.units import nbody_system
from amuse.units import units

from pylab import *

from datetime import date, timedelta

from amuse.support import data
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
        #self.setup_solar_system()
        #model t_0 should be the same as the 
        #loadstar t_0
        #OR 
        #pick t_0 from data for I.C.?
        self.model_t0 = date(1971,10,26)#CT 0.00
    

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

    def evolve(self, start_date, end_date, step):

        for x in arange(start_date, end_date, step):
            self.code.evolve_model((x-self.model_t0.toordinal())|units.day)
            if not self.code.get_indices_of_colliding_particles()['index_of_particle1'] == -1:
                pass
                
            #self.code.update_particles(self.planets)
            #self.planets.savepoint(x)

    def cleanup(self):
        del self.instance


def set_IC_at_date(data, particles, IC_date):

    for i, body in enumerate(data):
            
        particles[i].mass = body.mass
        particles[i].radius = body.radius
        r, v = body.get_vectors_at_date(IC_date)
        particles[i].position = units.AU([r[0],r[1],r[2]])
        particles[i].velocity = units.AUd([v[0],v[1],v[2]])

def approach(voyagerI, chkdistto):
    dR = voyagerI.position - chkdistto.position

    if (dot(dR,dR)**0.5) < 0.0080|units.AU:
        return True
    else:
        return False

def retreat(state):
    dist = (state['x']**2+state['y']**2+state['z']**2)**0.5
    if (dist|units.m).value_in(units.AU) > 0.0080:
        return True
    else:
        return False

if __name__ == '__main__':

    I = SolarSystemModel()

    planets_data = Horizons()

    #Notable events
    start_date = date(1971, 10, 26)
    voyagerI_launch_date = date(1977, 9, 7)
    stop_date = date(1989, 10, 26)

    I.model_t0 = start_date

    planets = data.Stars(10)
    voyagerI = data.Particle()

    #set I.C. using Horizons database
    set_IC_at_date([planets_data.Sun, planets_data.Mercury, planets_data.Venus, 
                    planets_data.Earth, planets_data.Moon, planets_data.Mars, 
                    planets_data.Jupiter, planets_data.Saturn, planets_data.Uranus,
                    planets_data.Neptune], 
                   planets,
                   start_date)

    set_IC_at_date([planets_data.VoyagerI],voyagerI.as_set(), voyagerI_launch_date)

    planets.synchronize_to(I.code.particles)
    
    planets_from_code_to_model = I.code.particles.new_channel_to(planets)
    
    #later we will add particles
    all_particles = planets.copy()

    P = planetarium.SolarSystemView((640,400))

    #Phase One #################################################################################
    for t in range(start_date.toordinal(), voyagerI_launch_date.toordinal(), 1):

        I.evolve(t,t+1,1)
        I.code.update_particles(all_particles)        
        
        planets_from_code_to_model.copy()

        P.renderamuse(all_particles)
        P.handle_events()
        if not P.go:
            P.go = True
            break
        
    #Phase Two #################################################################################
    all_particles.add_particle(voyagerI)
    all_particles.synchronize_to(I.code.particles)
    voyagerI_from_code_to_model = I.code.particles.new_channel_to(voyagerI.as_set())

    two_body_mode = False
    timestep = 1.0

    for t in arange(voyagerI_launch_date.toordinal(), stop_date.toordinal(), timestep):

        I.evolve(t,t+timestep,timestep)
        I.code.update_particles(all_particles)        
        
        planets_from_code_to_model.copy()
        
        if not two_body_mode:
            voyagerI_from_code_to_model.copy()
            for a_planet in planets:
                if approach(voyagerI, a_planet):
                    two_body_mode = True
                    print "close approach"
                    I.code.particles.remove_particle(voyagerI)
                    I.code.particles.remove_particle(a_planet) 
                    I.code.particles.synchronize_to(all_particles)
                    I.code.update_particles(all_particles)
                    #planets are heavy wtr spacecraft 
                    #so we put the planet back instead of the c.o.g. sys
                    all_particles.add_particle(a_planet)
                    all_particles.synchronize_to(I.code.particles)
                    a_planet_from_code_to_model =  I.code.particles.new_channel_to(a_planet.as_set())
                    I.code.update_particles(all_particles)
                    m = 1898.13e24#kg
                    #ra = 71492000.0#m
                    ra = 1.0#m

                    dr = voyagerI.position - a_planet.position
                    dv = voyagerI.velocity - a_planet.velocity
                    del voyagerI
                    tb = twobody.twobody()
                    tb.new_particle(m, 
                                    ra, 
                                    dr[0].value_in(units.m),
                                    dr[1].value_in(units.m),
                                    dr[2].value_in(units.m),
                                    dv[0].value_in(units.m/units.s),
                                    dv[1].value_in(units.m/units.s),
                                    dv[2].value_in(units.m/units.s))

                    days = 1
                    timestep = 0.05
                    break

        elif two_body_mode:
            
            tb.evolve(3600.0*24*1*days)
            days += timestep
            a_planet_from_code_to_model.copy()
            state, err = tb.get_state(0)
            
            if retreat(state):
                two_body_mode = False
                print "retreat"
                voyagerI = data.Particle()
                voyagerI_attrib = voyagerI.as_set()
                rel_vr = ([state['vx'],state['vy'],state['vz']]|units.m/units.s).as_quantity_in(units.AUd)
                rel_r = ([state['x'],state['y'],state['z']]|units.m).as_quantity_in(units.AU)
                voyagerI_attrib[0].velocity = a_planet.velocity + rel_vr
                voyagerI_attrib[0].position = a_planet.position + rel_r
                
                all_particles.add_particle(voyagerI)
                all_particles.synchronize_to(I.code.particles)
                voyagerI_from_code_to_model = I.code.particles.new_channel_to(voyagerI.as_set())
                #I.code.update_particles(all_particles)
                timestep = 1
                del tb
        
        P.renderamuse(all_particles)
        P.handle_events()
        if not P.go:
            P.go = True
            break
        

"""        
        if t==voyagerI_launch_date.toordinal():
            all_particles.add_particles(voyagerI)
            all_particles.synchronize_to(I.code.particles)
            voyager_from_code_to_model = I.code.particles.new_channel_to(voyagerI)


        if t>voyagerI_launch_date.toordinal():
            voyager_from_code_to_model.copy()
            #colli
            for number, chkdistto in enumerate(planets):
                dR = voyagerI.position - chkdistto.position
                if (dot(dR[0],dR[0])**0.5) < 0.1|units.AU:
                    #remove voyager and particle
                    I.code.remove_particles(voyagerI)                    
                    #add barycenter of pair
                    #evolve pair as twobody problem
"""
"""
    #add voyager to model
    I.voyagers = data.Stars(1)

    voyagerI = voyagers[0]

    #set I.C. for voyager 0
    voyagerI.mass = units.kg(500)
    voyagerI.radius = units.m(20)
    r, v = I.voyagerI.get_vectors_at_date(voyagerI_launch_date)
    voyagerI.position = units.AU([r[0],r[1],r[2]])
    voyagerI.velocity = units.AUd([v[0],v[1],v[2]])

    I.stars.add_particles(voyagers)
    I.stars.synchronize_to(I.code.particles)
    #I.instance.timestep = 0.01 

    #I.instance.setup_module()
    #I.instance.setup_particles(I.stars)

    two_body_mode = False

    for t in arange(voyagerI_launch_date.toordinal(),stop_date.toordinal(),1):
        r, v = voyagerI_data.get_vectors_at_date(date.fromordinal(t))
        DR = I.stars[10].position - units.AU([r[0],r[1],r[2]])
        yy.append((dot(DR,DR)**0.5).value_in(units.AU))
        I.evolve(t, t+1, 1)
        #Check if voyager comes close to any object:
        for i, s in enumerate(I.stars):
            if i in [6]:
                
                dr = (s.position - I.stars[10].position)
                dist = dot(dr,dr)**0.5

                if dist < two_body_mode_radius_ratio * s.radius:
                    if not two_body_mode:
                        print "two body integration with: "+bodies[i].name
                        print "close encounter, invoking two body integrator"
                        tb = twobody.twobody()
                        m = 1898.13e24#kg
                        #ra = 71492000.0#m
                        ra = 1.0#m
                        dv = (s.velocity-I.stars[10].velocity)

                        tb.new_particle(m, 
                                        ra, 
                                        dr[0].value_in(units.m),
                                        dr[1].value_in(units.m),
                                        dr[2].value_in(units.m),
                                        dv[0].value_in(units.m/units.s),
                                        dv[1].value_in(units.m/units.s),
                                        dv[2].value_in(units.m/units.s))
                        two_body_mode = True

                    tb.evolve(3600.0*24*1)
                    print ".",

                    state, err = tb.get_state(0)
                    I.stars[10].velocity[0] = ((state['vx']) | units.AUd) + I.stars[6].velocity[0]
                    I.stars[10].velocity[1] = ((state['vy']) | units.AUd) + I.stars[6].velocity[1]
                    I.stars[10].velocity[2] = ((state['vz']) | units.AUd) + I.stars[6].velocity[2]
                    I.stars[10].position[0] = ((state['x']) | units.AU) + I.stars[6].position[0]
                    I.stars[10].position[1] = ((state['y']) | units.AU) + I.stars[6].position[1]
                    I.stars[10].position[2] = ((state['z']) | units.AU) + I.stars[6].position[2]
                    
                    #I.stars.synchronize_to(I.instance.particles)
                        
                elif dist > two_body_mode_radius_ratio * s.radius:
                    if two_body_mode:
                        #assuming leaving close encounter area
                        print "Quitting two body mode"
                        two_body_mode = False
                        #I.stars.synchronize_to(I.instance.particles)

        I.stars.synchronize_to(I.instance.particles)
        P.renderamuse(I.stars)
        P.handle_events()
        if not P.go:
            P.go = True
            plot(yy,'o')
            show()
            s = raw_input()
            break

    S = []
"""
