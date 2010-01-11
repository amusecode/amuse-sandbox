import sys

from amuse.legacy.bhtree.interface import BHTreeInterface, BHTree
from amuse.legacy.hermite0.interface import HermiteInterface, Hermite
from amuse.legacy.phiGRAPE.interface import PhiGRAPEInterface, PhiGRAPE

from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
"""
JDCT 
   X     Y     Z
   VX    VY    VZ
   LT    RG    RR
Sun-mercury
$$SOE
2454797.500000000 = A.D. 2008-Nov-27 00:00:00.0000 (CT)
  -1.812507519383936E-01 -4.238556570722329E-01 -1.858336536398257E-02
   2.028905659997442E-02 -9.482619060967475E-03 -2.636707283074494E-03
   2.664579022492631E-03  4.613575561087286E-01  8.471820990251489E-04

Sun-Venus
$$SOE
2454797.500000000 = A.D. 2008-Nov-27 00:00:00.0000 (CT)
   7.174859394658725E-01 -9.118094489213757E-02 -4.286369239375957E-02
   2.569149540621095E-03  1.995467986481682E-02  1.246915402703626E-04
   4.184510840750728E-03  7.245255924867496E-01  2.553031535710333E-05

Earth-sun
2454797.500000000 = A.D. 2008-Nov-27 00:00:00.0000 (CT)
   4.152751345538410E-01  8.988236789078493E-01 -4.821560120168533E-05
  -1.588315644843389E-02  7.210443860976745E-03 -1.350251461569625E-07
   5.718455717278124E-03  9.901199146915970E-01 -1.161094180263466E-04
Sun-moon
2454797.500000000 = A.D. 2008-Nov-27 00:00:00.0000 (CT
   4.138191074397691E-01  8.965602570292573E-01 -2.762446418149536E-04
  -1.541527312550390E-02  6.894586206029982E-03  1.223837010915995E-05
   5.703062252841441E-03  9.874546189459633E-01 -2.002380278932658E-04
Sun-Mars
2454797.500000000 = A.D. 2008-Nov-27 00:00:00.0000 (CT)
  -5.350147323170708E-01 -1.400272441929516E+00 -1.637545552747233E-02
   1.360625065710822E-02 -3.765876077406818E-03 -4.130340644254660E-04
   8.658023706817600E-03  1.499090334491985E+00 -1.333827852316589E-03
Sun-Jupiter
2454797.500000000 = A.D. 2008-Nov-27 00:00:00.0000 (CT)
   2.503092399973117E+00 -4.468134118102924E+00 -3.752173268244928E-02
   6.490840561446090E-03  4.046895067472646E-03 -1.620422227298534E-04
   2.958007250792838E-02  5.121630789170776E+00 -3.570769184429194E-04

Sun-Saturn
2454797.500000000 = A.D. 2008-Nov-27 00:00:00.0000 (CT)
  -9.023156820924056E+00  2.468810475231705E+00  3.161126539154331E-01
  -1.769097295887704E-03 -5.393257979873611E-03  1.639859191780030E-04
   5.405968809328438E-02  9.360144837958844E+00  2.884302120086998E-04

"""
class LoadStar(object):

    def __init__(self, planet):
        f = open(planet+'.txt','r')
        M = f.readlines()
        f.close()

        radiusline = M[4]
        radius = float(radiusline.split('=')[1].split('(')[0])

        massline =  M[5]
        Myunit = massline.split('(')[1].split(')')[0].split('^')
        Unitbase = float(Myunit[0])
        Expo = float(Myunit[1].split('kg')[0])
        mass = Unitbase**Expo

        days = []

        for i, s in enumerate(M):
            if 'A.D.' in s:
                days.append(i)

        days.pop(0)#start date
        days.pop(0)#stop date

        
        r = []
        v = []
        for i in days:
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
            
        self.r = r
        self.v = v
        self.mass = mass
        self.radius = radius

class SolarSystemModel(object):

    def setup_solar_system(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        
        stars = core.Stars(8)
        
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = units.AU(array((0.0,0.0,0.0)))
        sun.velocity = units.AUd(array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)

   
        mercury = stars[1]
        mercury.mass = units.kg( 3.302e23)
        mercury.radius = units.km(2440) 
        mercury.position = units.AU(array((-1.812507519383936E-01, -4.238556570722329E-01, -1.858336536398257E-02)))
        mercury.velocity = units.AUd(array((2.028905659997442E-02, -9.482619060967475E-03, -2.636707283074494E-03)))

        venus = stars[2]
        venus.mass = units.kg(48.685e23)
        venus.radius = units.km(6051.8) 

        venus.position = units.AU(array((7.174859394658725E-01, -9.118094489213757E-02, -4.286369239375957E-02)))
        venus.velocity = units.AUd(array((2.569149540621095E-03,  1.995467986481682E-02,  1.246915402703626E-04)))

        earth = stars[3]
        earth.mass = units.kg(5.9736e24)

        earth.radius = units.km(6371) 
        earth.position = units.AU(array((   4.152751345538410E-01,  8.988236789078493E-01, -4.821560120168533E-05)))
        earth.velocity = units.AUd(array((  -1.588315644843389E-02,  7.210443860976745E-03, -1.350251461569625E-07)))
      
        moon = stars[4]
        moon.mass = units.kg(7.35e22)
        moon.radius = units.km(1738) 
        moon.position = units.AU(array((   4.138191074397691E-01,  8.965602570292573E-01, -2.762446418149536E-04)))
        moon.velocity = units.AUd(array((  -1.541527312550390E-02,  6.894586206029982E-03,  1.223837010915995E-05))) 
        
        mars = stars[5]
        mars.mass = units.kg(6.4185e23)
        mars.radius = units.km(3389.9)
        mars.position = units.AU(array((  -5.350147323170708E-01, -1.400272441929516E+00, -1.637545552747233E-02)))
        mars.velocity = units.AUd(array((  1.360625065710822E-02, -3.765876077406818E-03, -4.130340644254660E-04))) 
        
        jupiter = stars[6]
        jupiter.mass = units.kg(1898.13e24)
        jupiter.radius = units.km(71492)
        jupiter.position = units.AU(array((   2.503092399973117E+00, -4.468134118102924E+00, -3.752173268244928E-02)))
        jupiter.velocity = units.AUd(array((   6.490840561446090E-03,  4.046895067472646E-03, -1.620422227298534E-04))) 

        saturn = stars[7]
        saturn.mass = units.kg(5.68319e26 )
        saturn.radius = units.km(58232)
        saturn.position = units.AU(array((   -9.023156820924056E+00,  2.468810475231705E+00,  3.161126539154331E-01)))
        saturn.velocity = units.AUd(array(( -1.769097295887704E-03, -5.393257979873611E-03,  1.639859191780030E-04  ))) 

        dist = earth.position - moon.position 
        print "distance vector"+str(dist)

        print "distance %s\n" %  (dot(dist,dist)**0.5).as_quantity_in(units.m)#(dist*dist).sum()**0.5 

        velo = moon.velocity-earth.velocity
        print "orb velocity %s\n" % (dot(velo,velo)**0.5).as_quantity_in(units.m/units.s)

        return stars

    def evolve(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        #instance = PhiGRAPE(convert_nbody)
        instance = BHTree(convert_nbody)
        #instance.dt_dia = 1
        instance.parameters.epsilon_squared = 0.00000001 | units.AU**2
        #instance.set_eta(0.0001,0.0002)
        instance.dt_dia = 5000
        instance.timestep = 0.0001 
        #instance.use_self_gravity = 0
        instance.setup_module()

        stars = self.setup_solar_system()

        instance.setup_particles(stars)
        #instance.initialize_particles(0.0)

        s=raw_input()
          
        ion()
        for x in arange(1,2.0*365,10):
            instance.evolve_model(x|units.day)
            tijd, error =  instance.get_time()
            nbodytijd =  tijd|nbody_system.time
            print convert_nbody.to_si(nbodytijd).as_quantity_in(units.day)
          
            energie, error = instance.get_kinetic_energy()

            print energie
            instance.update_particles(stars)
            stars.savepoint()


        fig = figure()
        #plotje = fig.add_subplot(1,1,1, aspect='equal')
        
        ax = Axes3D(fig)

        for i in range(len(stars)):
            x_points = stars[i].get_timeline_of_attribute("x")
            y_points = stars[i].get_timeline_of_attribute("y")
            
            x_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), x_points)
            y_points_in_AU = map(lambda (t,x) : x.value_in(units.AU), y_points)
            
            ax.plot(x_points_in_AU,y_points_in_AU, linewidth = 0.5)#, color = "k", marker = '+')
            
            ax.set_xlim3d(-20, 20)
            ax.set_ylim3d(-20, 20)
            ax.set_zlim3d(-20, 20)
        
        
        show()
        s = raw_input()
        fig.savefig("solsys.svg")    
        instance.cleanup_module()
        del instance

if __name__ == '__main__':
  Mercury = LoadStar('mercury')
  print Mercury.r[0][0]
  print Mercury.mass
  #I = SolarSystemModel()
  #I.evolve()

























