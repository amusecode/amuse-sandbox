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

from amuse import datamodel
class SolarSystemModel(object):

    def __init__(self):
        self.load_solar_system()
        self.load_integrator_and_units()
        self.stars = []
        #self.setup_solar_system()
        #model t_0 should be the same as the 
        #loadstar t_0
        #OR 
        #pick t_0 from data for I.C.?
        self.model_t0 = date(1971,10,26)#CT 0.00
    
    def load_solar_system(self):

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
        self.VoyagerI.mass = units.kg(1000)
        self.VoyagerI.radius = units.m(20)

        self.VoyagerII = LoadStar('voyagerII')
        self.VoyagerII.mass = units.kg(1000)
        self.VoyagerII.radius = units.m(20)

    def load_integrator_and_units(self):
        self.convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        #self.instance = PhiGRAPE(self.convert_nbody)
        #self.instance = BHTree(self.convert_nbody)
        self.instance = Hermite(self.convert_nbody)
        self.instance.parameters.epsilon_squared = 0.000000001 | units.AU**2
        #self.instance.set_eta(0.0001,0.0002)
        self.instance.dt_dia = 10
        self.instance.timestep = 0.0002 
        #instance.use_self_gravity = 0
        self.instance.setup_module()

    def evolve(self, start_date, end_date, step):

        stars = self.stars
        #instance.initialize_particles(0.0)

        for x in arange(start_date, end_date, step):
            self.instance.evolve_model((x-self.model_t0.toordinal())|units.day)
            tijd, error = self.instance.get_time()
            nbodytijd =  tijd | nbody_system.time
            #print self.convert_nbody.to_si(nbodytijd).as_quantity_in(units.day)
            energie, error = self.instance.get_kinetic_energy()
            #print energie
            """
            if len(stars)>10:
                #assumde added voyagers
                vdist = stars[6].position-stars[10].position
                dist =  (dot(vdist,vdist)**0.5).as_quantity_in(units.km)
                if dist<1e6|units.km:
                    print "close to Jupiter"
            """        
            self.instance.update_particles(stars)

            stars.savepoint(x)

        #self.instance.cleanup_module()
        #del instance we do that somewhere else...

    def distance(self, x1,y1,z1, x2,y2,z2):
        dx = x2-x1
        dy = y2-y1
        dz = z2-z1
        return (dx**2+dy**2+dz**2)**0.5

    def cleanup(self):
        del self.instance

if __name__ == '__main__':

    I = SolarSystemModel()
    
    #make bodies in mem
    I.stars = datamodel.Stars(10)

    two_body_mode_radius = 300

    bodies = [I.Sun, I.Mercury, I.Venus, I.Earth, I.Moon, I.Mars, 
              I.Jupiter, I.Saturn, I. Uranus,I.Neptune]

    #set IC for bodies in mem according to data set and date
    start_date = date(1971, 10, 26)
    voyagerI_launch_date = date(1977, 9, 7)

    for i, body in enumerate(bodies):
            
        I.stars[i].mass = body.mass
        I.stars[i].radius = body.radius
        r, v = body.get_vectors_at_date(start_date)
        I.stars[i].position = units.AU(array((r[0],r[1],r[2])))
        I.stars[i].velocity = units.AUd(array((v[0],v[1],v[2])))

    #I.instance.setup_particles(I.stars)
    I.stars.synchronize_to(I.instance.particles)
        
    I.model_t0 = start_date#date(1971,10,26)    

    P = planetarium.SolarSystemView((640,400))

    for t in range(start_date.toordinal(), voyagerI_launch_date.toordinal(),1):
        I.evolve(t,t+1,1)
        P.renderamuse(I.stars)
        P.handle_events()
        if not P.go:
            P.go = True
            break

    voyagers = datamodel.Stars(1)
    voyagerI = voyagers[0]
    VoyagerI = LoadStar('voyagerI')
    voyagerI.mass = units.kg(500)
    voyagerI.radius = units.m(20)
    r, v = VoyagerI.get_vectors_at_date(voyagerI_launch_date)
    voyagerI.position = units.AU(array((r[0],r[1],r[2])))
    voyagerI.velocity = units.AUd(array((v[0],v[1],v[2])))
    I.stars.add_particles(voyagers)
    I.stars.synchronize_to(I.instance.particles)
    I.instance.timestep = 0.002 
    #I.instance.setup_module()
    #I.instance.setup_particles(I.stars)

    two_body_mode = False

    for t in arange(voyagerI_launch_date.toordinal(),date(1987,9,7).toordinal(),0.2):
        I.evolve(t, t+1, 0.2)#till voyager launch + some days
        #Check if voyager comes close to any object:
        for i, s in enumerate(I.stars):
            if i in [6]:
                
                dx = (s.position[0]-I.stars[10].position[0]).value_in(units.AU)
                dy = (s.position[0]-I.stars[10].position[0]).value_in(units.AU)
                dz = (s.position[0]-I.stars[10].position[0]).value_in(units.AU)
                dist =  (dx**2+dy**2+dz**2)**0.5
                if dist<two_body_mode_radius * s.radius.value_in(units.AU):
                    if not two_body_mode:
                        print "two body integration with: "+bodies[i].name
                        print "close encounter, invoking two body integrator"
                        tb = twobody.twobody()
                        m = 1898.13e24#kg
                        ra = 71492000.0#m
                        x = (dx|units.AU).value_in(units.m)
                        y = (dy|units.AU).value_in(units.m)
                        z = (dz|units.AU).value_in(units.m)
                        dvx = (s.velocity[0]-I.stars[10].velocity[0]).value_in(units.AUd)
                        dvy = (s.velocity[1]-I.stars[10].velocity[1]).value_in(units.AUd)
                        dvz = (s.velocity[2]-I.stars[10].velocity[2]).value_in(units.AUd)
                        vx= (dvx|units.AUd).value_in(units.m/units.s)
                        vy= (dvy|units.AUd).value_in(units.m/units.s)
                        vz= (dvz|units.AUd).value_in(units.m/units.s)
                        tb.new_particle(m, ra, x,y,z,vx,vy,vz)
                        two_body_mode = True
                    #tb.evolve((t-voyagerI_launch_date.toordinal())*3600.0*24)
                    tb.evolve(3600.0*24*0.2)
                    state, err = tb.get_state(0)
                    I.stars[10].velocity[0] = ((state['vx']) | units.AUd) + I.stars[6].velocity[0]
                    I.stars[10].velocity[1] = ((state['vy']) | units.AUd) + I.stars[6].velocity[1]
                    I.stars[10].velocity[2] = ((state['vz']) | units.AUd) + I.stars[6].velocity[2]
                    I.stars[10].position[0] = ((state['x']) | units.AU) + I.stars[6].position[0]
                    I.stars[10].position[1] = ((state['y']) | units.AU) + I.stars[6].position[1]
                    I.stars[10].position[2] = ((state['z']) | units.AU) + I.stars[6].position[2]
                    I.stars.synchronize_to(I.instance.particles)
                        
                elif dist>two_body_mode_radius * s.radius.value_in(units.AU):
                    if two_body_mode:
                        #assuming leaving close encounter area
                        print "Quitting two body mode"
                        two_body_mode = False
                        I.stars.synchronize_to(I.instance.particles)


        P.renderamuse(I.stars)
        P.handle_events()
        if not P.go:
            P.go = True
            break

    S = []
