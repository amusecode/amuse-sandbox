import sys

from amuse.legacy.bhtree.interface import BHTreeInterface, BHTree
from amuse.legacy.hermite0.interface import HermiteInterface, Hermite
from amuse.legacy.phiGRAPE.interface import PhiGRAPEInterface, PhiGRAPE
import amuse.legacy.twobody.twobody as twobody

from amuse.units import nbody_system
from amuse.units import units
from mpl_toolkits.mplot3d import Axes3D
from pylab import *

from datetime import date, timedelta

from amuse import datamodel
class LoadStar(object):

    def __init__(self, planet):
        f = open(planet+'.txt','r')
        M = f.readlines()
        f.close()
        """
        radiusline = M[4]
        radius = float(radiusline.split('=')[1].split('(')[0])

        massline =  M[5]
        Myunit = massline.split('(')[1].split(')')[0].split('^')
        Unitbase = float(Myunit[0])
        Expo = float(Myunit[1].split('kg')[0])
        mass = Unitbase**Expo
        """
        days = []

        for i, s in enumerate(M):
            if 'A.D.' in s:
                days.append(i)

        days.pop(0)#start date
        days.pop(0)#stop date

        
        r = []
        v = []
        julian_day = []
        
        for i in days:
            
            julian_day.append((float(M[i].split('=')[0])-1721424.5))

            #http://ssd.jpl.nasa.gov/?horizons_doc#time
            rs = M[i+1].split(' ')
            vs = M[i+2].split(' ')
            foor = []
            foov = []
            for j, n in enumerate(rs):
                if not n == '':
                    foor.append(float(n))
            for j, n in enumerate(vs):
                if not n == '':
                    foov.append(float(n))
            
                    
            r.append([foor[0], foor[1], foor[2]])
            v.append([foov[0], foov[1], foov[2]])
            
        self.ordinal = julian_day
        self.r = r
        self.v = v
        #not implemented yet:
        self.mass = 0.0
        self.radius = 0.0
        del M

    def get_vectors_at_date(self, at_date):
        indices = [i for i, j in enumerate(self.ordinal) if j==at_date.toordinal()]
        return self.r[indices[0]], self.v[indices[0]]

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
            tijd, error =  self.instance.get_time()
            nbodytijd =  tijd | nbody_system.time
            #print self.convert_nbody.to_si(nbodytijd).as_quantity_in(units.day)
            energie, error = self.instance.get_kinetic_energy()
            #print energie

            if len(stars)>10:
                #assumde added voyagers
                vdist = stars[6].position-stars[10].position
                dist =  (dot(vdist,vdist)**0.5).as_quantity_in(units.km)
                if dist<1e6|units.km:
                    print "close to Jupiter"

            self.instance.update_particles(stars)
            stars.savepoint()

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
    I.evolve(date(1971,10,26),date(1978,1,3),10)#till voyager launch + some days
    
    voyagers = datamodel.Stars(1)

    voyagerI = voyagers[0]#I.stars.new_particle()
    #voyagerII = voyagers[1]# I.stars.new_particle()
 
    VoyagerI = LoadStar('voyagerI')
    voyagerI.mass = units.kg(1000)
    voyagerI.radius = units.m(20)
    r, v = VoyagerI.get_vectors_at_date(date(1978,1,3))

    Jx = I.stars[6].position.x.value_in(units.m)
    Jy = I.stars[6].position.y.value_in(units.m)
    Jz = I.stars[6].position.z.value_in(units.m)
    Jvx = I.stars[6].velocity.x.value_in(units.m/units.s)
    Jvy = I.stars[6].velocity.y.value_in(units.m/units.s)
    Jvz = I.stars[6].velocity.z.value_in(units.m/units.s)

    
    s = raw_input()
    voyagerI.position = units.AU(array((r[0],r[1],r[2])))
    voyagerI.velocity = units.AUd(array((v[0],v[1],v[2])))
    """
    VoyagerII = LoadStar('voyagerII')    
    voyagerII.mass = units.kg(1000)
    voyagerII.radius = units.m(20)
    r, v = VoyagerII.get_vectors_at_date(date(1977,7,3))
    voyagerII.position = units.AU(array((r[0],r[1],r[2])))
    voyagerII.velocity = units.AUd(array((v[0],v[1],v[2])))
    """


    I.stars.add_particles(voyagers)
    
    I.stars.synchronize_to(I.instance.particles)

    I.instance.timestep = 0.001 
    
    I.instance.setup_module()
    I.instance.setup_particles(I.stars)
    I.evolve(date(1978,1,3),date(1980,9,5),1)
    
    
    fig = figure()
    ax = Axes3D(fig)
    ion()

    for i in range(len(I.stars)):
        x_points = I.stars[i].get_timeline_of_attribute("x")
        y_points = I.stars[i].get_timeline_of_attribute("y")
        z_points = I.stars[i].get_timeline_of_attribute("z")
        x_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), x_points)
        y_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), y_points)
        z_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), z_points)
        ax.plot(x_points_in_AU,y_points_in_AU, z_points_in_AU, linewidth = 0.5)#, color = "k", marker = '+')

    nb=twobody.twobody()
    m = 1898.13e24#kg
    ra = 71492000.0#m
    x = (r[0]|units.AU).value_in(units.m)-Jx
    y = (r[1]|units.AU).value_in(units.m)-Jy
    z = (r[2]|units.AU).value_in(units.m)-Jz
    vx= (v[0]|units.AUd).value_in(units.m/units.s)-Jvx
    vy= (v[1]|units.AUd).value_in(units.m/units.s)-Jvy
    vz= (v[2]|units.AUd).value_in(units.m/units.s)-Jvz
    print m, ra, x, y, z, vx, vy, vz
    nb.new_particle(m, ra, x, y, z, vx, vy, vz)
    g
    #examples of how to update the data of one particle or all of 'em

    #I.stars[6].position = [x,y,z]|unit.AU
    
    #code.particles[6].postition = ...
 
    #I.stars.copy_values_of_state_attributes_to(code.particles)
    #I.stars.copy_values_of_state_attributes_to(code.particles)

    x_points = I.stars[6].get_timeline_of_attribute("x")
    y_points = I.stars[6].get_timeline_of_attribute("y")
    z_points = I.stars[6].get_timeline_of_attribute("z")
    x_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), x_points)
    y_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), y_points)
    z_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), z_points)
    
    print len(x_points_in_AU)

    for i,j  in enumerate(range(0,3600.0*24*180.13, 3600.0*24)):
        print i
        nb.evolve(j)
        state,err=nb.get_state(0)
        ax.scatter( 
                   [((state['x']) | units.m).value_in(units.AU)+x_points_in_AU[i]], 
                   [((state['y']) | units.m).value_in(units.AU)+y_points_in_AU[i]],
                   [((state['z']) | units.m).value_in(units.AU)+z_points_in_AU[i]],
                 marker = '+')
    
    #ax.plot(x1,y1,z1)
    ax.set_xlim3d(-10, 10)
    ax.set_ylim3d(-10, 10)
    ax.set_zlim3d(-10, 10) 
    
    show()

"""    
    x1 = [i[0] for i in VoyagerI.r]
    y1 = [i[1] for i in VoyagerI.r]
    z1 = [i[2] for i in VoyagerI.r]

    def setup_solar_system(self):
        
        print "Loading ephemeris of solar system objects from Horizons"

        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        
        stars = datamodel.Stars(10)
        
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = units.AU(array((0.0,0.0,0.0)))
        sun.velocity = units.AUd(array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)
  
        mercury = stars[1]
        mercury.mass = units.kg( 3.302e23)
        mercury.radius = units.km(2440) 
        mercury.position = units.AU(array(self.Mercury.r[0]))
        mercury.velocity = units.AUd(array(self.Mercury.v[0]))

        venus = stars[2]
        venus.mass = units.kg(48.685e23)
        venus.radius = units.km(6051.8) 

        venus.position = units.AU(array((self.Venus.r[0])))
        venus.velocity = units.AUd(array((self.Venus.v[0])))

        earth = stars[3]
        earth.mass = units.kg(5.9736e24)

        earth.radius = units.km(6371) 
        earth.position = units.AU(array((self.Earth.r[0])))
        earth.velocity = units.AUd(array((self.Earth.v[0])))
      
        moon = stars[4]
        moon.mass = units.kg(7.35e22)
        moon.radius = units.km(1738) 
        moon.position = units.AU(array((self.Moon.r[0])))
        moon.velocity = units.AUd(array((self.Moon.v[0]))) 
        
        mars = stars[5]
        mars.mass = units.kg(6.4185e23)
        mars.radius = units.km(3389.9)
        mars.position = units.AU(array((self.Mars.r[0])))
        mars.velocity = units.AUd(array((self.Mars.v[0]))) 
        
        jupiter = stars[6]
        jupiter.mass = units.kg(1898.13e24)
        jupiter.radius = units.km(71492)
        jupiter.position = units.AU(array((self.Jupiter.r[0])))
        jupiter.velocity = units.AUd(array((self.Jupiter.v[0]))) 

        saturn = stars[7]
        saturn.mass = units.kg(5.68319e26 )
        saturn.radius = units.km(58232)
        saturn.position = units.AU(array((self.Saturn.r[0])))
        saturn.velocity = units.AUd(array((self.Saturn.v[0]))) 

        uranus = stars[8]
        uranus.mass = units.kg(86.8103e24 )
        uranus.radius = units.km(26000)
        uranus.position = units.AU(array((self.Uranus.r[0])))
        uranus.velocity = units.AUd(array((self.Uranus.v[0]))) 

        neputne = stars[9]
        neputne.mass = units.kg(102.41e24 )
        neputne.radius = units.km(25000)
        neputne.position = units.AU(array((self.Neptune.r[0])))
        neputne.velocity = units.AUd(array((self.Neptune.v[0]))) 
        
        dist = earth.position - moon.position 

        #print "distance vector"+str(dist)
        #print "distance %s\n" %  (dot(dist,dist)**0.5).as_quantity_in(units.m)#(dist*dist).sum()**0.5 
        #velo = moon.velocity-earth.velocity
        #print "orb velocity %s\n" % (dot(velo,velo)**0.5).as_quantity_in(units.m/units.s)

        self.stars = stars
    def plot_sim(self, range1, range2):
        
        fig = figure()
        ax = Axes3D(fig)
  
        ion()

        stars = self.stars

        for i in range(len(stars)):
            x_points = stars[i].get_timeline_of_attribute("x")
            y_points = stars[i].get_timeline_of_attribute("y")
            z_points = stars[i].get_timeline_of_attribute("z")
            x_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), x_points)
            y_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), y_points)
            z_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), z_points)
            ax.plot(x_points_in_AU,y_points_in_AU, z_points_in_AU, linewidth = 0.5)#, color = "k", marker = '+')
            
        for i in range2:
            #ax.plot([self.Mercury.r[i][0]],[self.Mercury.r[i][1]],[self.Mercury.r[i][2]],'o')
            ax.plot([self.Earth.r[i][0]],[self.Earth.r[i][1]],[self.Earth.r[i][2]],'o')
            ax.plot([self.Jupiter.r[i][0]],[self.Jupiter.r[i][1]],[self.Jupiter.r[i][2]],'o')
            ax.plot([self.Saturn.r[i][0]],[self.Saturn.r[i][1]],[self.Saturn.r[i][2]],'o')
            ax.plot([self.Uranus.r[i][0]],[self.Uranus.r[i][1]],[self.Uranus.r[i][2]],'o')
            ax.plot([self.Neptune.r[i][0]],[self.Neptune.r[i][1]],[self.Neptune.r[i][2]],'o')

        print self.Mercury.r[-1][0]/x_points_in_AU[-1]
        print self.Mercury.r[-1][1]/y_points_in_AU[-1]
        print self.Mercury.r[-1][2]/z_points_in_AU[-1]
        
        ax.set_xlim3d(-10, 10)
        ax.set_ylim3d(-10, 10)
        ax.set_zlim3d(-10, 10) 
        #show()
        fig.savefig("solsys.svg")    
        return ax


"""
