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
        self.instance = BHTree(self.convert_nbody)
        #self.instance = Hermite(self.convert_nbody)
        #self.instance.dt_dia = 1
        self.instance.parameters.epsilon_squared = 0.000000001 | units.AU**2
        #instance.set_eta(0.0001,0.0002)
        self.instance.dt_dia = 5000
        self.instance.timestep = 0.001 
        #instance.use_self_gravity = 0
        self.instance.setup_module()

    def evolve(self, start_date, end_date, step):

        stars = self.stars
        #instance.initialize_particles(0.0)

        for x in arange(start_date.toordinal(), end_date.toordinal(), step):
            self.instance.evolve_model((x-self.model_t0.toordinal())|units.day)
            tijd, error = self.instance.get_time()
            nbodytijd =  tijd | nbody_system.time
            #print self.convert_nbody.to_si(nbodytijd).as_quantity_in(units.day)
            energie, error = self.instance.get_kinetic_energy()
            #print energie
            print stars[0].position
            if len(stars)>10:
                #assumde added voyagers
                vdist = stars[6].position-stars[10].position
                dist =  (dot(vdist,vdist)**0.5).as_quantity_in(units.km)
                if dist<1e6|units.km:
                    print "close to Jupiter"

            self.instance.update_particles(stars)

            stars.savepoint(x)

        #self.instance.cleanup_module()
        #del instance we do that somewhere else...

    def cleanup(self):
        del self.instance


if __name__ == '__main__':
    

    I = SolarSystemModel()
    
    #make bodies in mem
    I.stars = datamodel.Stars(10)

    bodies = [I.Sun, I.Mercury, I.Venus, I.Earth, I.Moon, I.Mars, 
              I.Jupiter, I.Saturn, I. Uranus,I.Neptune]

    #set IC for bodies in mem according to data set and date
    for i, body in enumerate(bodies):
            
        I.stars[i].mass = body.mass
        I.stars[i].radius = body.radius
        r, v = body.get_vectors_at_date(date(1971, 10, 26))
        I.stars[i].position = units.AU(array((r[0],r[1],r[2])))
        I.stars[i].velocity = units.AUd(array((v[0],v[1],v[2])))

    I.instance.setup_particles(I.stars)
    
    I.model_t0 = date(1971,10,26)    
    I.evolve(date(1971,10,26),date(1977,7,1)+timedelta(100),1)#till voyager launch + some days
    
    P = planetarium.SolarSystemView((800,600))
    S = []
    
    for i in range(len(I.stars)):
        S.append(NewStar())
        x_points = I.stars[i].get_timeline_of_attribute("x")
        y_points = I.stars[i].get_timeline_of_attribute("y")
        z_points = I.stars[i].get_timeline_of_attribute("z")

        x_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), x_points)
        y_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), y_points)
        z_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), z_points)
        #make date array

        S[-1].ordinal = map(lambda (t,x) : t, x_points)
        print S[-1].ordinal
        #for some reason last timestamp is None !!!!
        S[-1].max = S[-1].ordinal[-2]
        S[-1].r = array([x_points_in_AU, y_points_in_AU, z_points_in_AU]).transpose()
        S[-1].v = array([x_points_in_AU, y_points_in_AU, z_points_in_AU]).transpose()
        P.add_planet(S[-1])
        
    while P.go:
        
        P.handle_events()
        P.animate()



