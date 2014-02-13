import numpy

#from amuse.community.test_kepler.interface import Kepler
from interface import Kepler


from amuse.units import nbody_system
from amuse.units import units,constants

from amuse.ic.plummer import new_plummer_model

from amuse.datamodel import Particle

from matplotlib import pyplot

import time

from amuse.ext.orbital_elements import orbital_elements_from_binary

def elements(starmass,x,y,z,vx,vy,vz,G=constants.G):
    mu=G*starmass
    r=(x**2+y**2+z**2)**0.5
    v2=(vx**2+vy**2+vz**2)
    
    e=v2/2-mu/r
    
    a=-mu/2/e
    
    hx=y*vz-z*vy
    hy=z*vx-x*vz
    hz=x*vy-y*vx

    rdotv=x*vx+y*vy+z*vz

    ex=v2*x/mu-rdotv*vx/mu-x/r
    ey=v2*y/mu-rdotv*vy/mu-y/r
    ez=v2*z/mu-rdotv*vz/mu-z/r

    h2=(hx**2+hy**2+hz**2)    
    
    eps=(1-h2/mu/a)**0.5
    
    return a,eps

def test_kepler( N,tend=1.| units.yr,method=0):

  numpy.random.seed(1234567)
  
  conv=nbody_system.nbody_to_si(2.| units.MSun, 5.|units.AU)
  
  comets=new_plummer_model(N,conv)
  
  sun=Particle(mass=1.|units.MSun)
  
  sun.position=[0,0,0]|units.AU
  sun.velocity=[0,0,0]|units.kms
  
  comets.mass*=0.
  
  code=Kepler(conv,redirection="none")
  
  code.set_method(method)
  
  code.central_particle.add_particle(sun)
  code.orbiters.add_particles(comets)
  
  a0,eps0=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz)
  
#  pyplot.plot(a0.value_in(units.AU),eps0,"ro")
  
#  print code.orbiters.x[0]
  print orbital_elements_from_binary(code.particles[0:2],constants.G)
  
  t1=time.time()
  code.evolve_model(tend)
  t2=time.time()
  
  print orbital_elements_from_binary(code.particles[0:2],constants.G)
#  print code.orbiters.x[0]
  
  
  
  a,eps=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz)
  
#  pyplot.plot(a.value_in(units.AU),eps,"g+")
  
  da=(a-a0)/a0
  deps=(eps-eps0)/eps0
  
  print "max da,deps:",da.max(), deps.max()
  
  print "time:",t2-t1

#  pyplot.show()
  
  return t2-t1,da.max(),deps.max()
  
  
if __name__=="__main__":
  for method in [0,1]:
    test_kepler(N=10000,tend=1.e9| units.yr,method=method)  

#  for method in [0,1]:
#    test_kepler(N=100000,tend=1000.| units.yr,method=method)  


#  for method in [0,1]:
#    test_kepler(N=100,tend=1e6| units.yr,method=method)  
  
