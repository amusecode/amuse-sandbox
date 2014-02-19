import numpy

#from amuse.community.test_kepler.interface import Kepler
from interface import Kepler


from amuse.units import nbody_system
from amuse.units import units,constants

from amuse.ic.plummer import new_plummer_model

from amuse.datamodel import Particle

from matplotlib import pyplot

import time

from amuse.ext.orbital_elements import orbital_elements_from_binary,new_binary_from_orbital_elements

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
  
  
#  print code.orbiters.x[0]
  print orbital_elements_from_binary(code.particles[0:2],constants.G)
  
  t1=time.time()
  code.evolve_model(tend)
  t2=time.time()
  
  print orbital_elements_from_binary(code.particles[0:2],constants.G)
#  print code.orbiters.x[0]
  
  
  
  a,eps=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz)
  
  da=abs((a-a0)/a0)
  deps=abs(eps-eps0)/eps0

  dev=numpy.where(da > 0.00001)[0]
  print len(dev)

  print a0[dev].value_in(units.AU)
  print eps0[dev]

  pyplot.plot(a0[dev].value_in(units.AU),eps0[dev],"ro")
  pyplot.plot(a[dev].value_in(units.AU),eps[dev],"g+")
  
  
  print "max da,deps:",da.max(), deps.max()
  
  print "time:",t2-t1

  pyplot.show()
  
  return t2-t1,da.max(),deps.max()
  
def test_kepler_almost_parabolic( tend=1,method=0):
  code=Kepler(redirection="none")
  
  code.set_method(method)
  
  mass1=1.| nbody_system.mass
  mass2=0| nbody_system.mass
  semimajor_axis=1.|nbody_system.length
  eccentricity=0.999999
  p=2*numpy.pi*(semimajor_axis**3/nbody_system.G/mass1)**0.5
  tend=tend*p
  print tend
  parts=new_binary_from_orbital_elements(
          mass1,
          mass2,
          semimajor_axis, 
          eccentricity = eccentricity,
          true_anomaly = 0.0102121
      )

  code.central_particle.add_particle(parts[0])
  code.orbiters.add_particle(parts[1])

  a0,eps0=elements(mass1,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)
  
  
  print orbital_elements_from_binary(code.particles[0:2])
  
  t1=time.time()
  code.evolve_model(tend)
  t2=time.time()

  print "time:",t2-t1

  
  print orbital_elements_from_binary(code.particles[0:2])
  
  print code.orbiters.position

  
  a,eps=elements(mass1,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)
  
  da=abs((a-a0)/a0)
  deps=abs(eps-eps0)/eps0

  print da,deps

def test_kepler_parabolic( tend=1,method=0):
  code=Kepler(redirection="none")
  
  code.set_method(method)
  
  sun=Particle()
  sun.mass=1. | nbody_system.mass
  sun.x=0. | nbody_system.length
  sun.y=0. | nbody_system.length
  sun.z=0. | nbody_system.length
  sun.vx=0. | nbody_system.speed
  sun.vy=0. | nbody_system.speed
  sun.vz=0. | nbody_system.speed

  comet=Particle()
  comet.mass= 0 | nbody_system.mass
  comet.x=1. | nbody_system.length
  comet.y=0. | nbody_system.length
  comet.z=0. | nbody_system.length
  comet.vx=0. | nbody_system.speed
  comet.vy=1.0000000001*(2*nbody_system.G*sun.mass/comet.x)**0.5 
  comet.vz=0. | nbody_system.speed
  
  tend=tend | nbody_system.time
  print tend

  code.central_particle.add_particle(sun)
  code.orbiters.add_particle(comet)

  a0,eps0=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)
  
  print orbital_elements_from_binary(code.particles[0:2])
  
  t1=time.time()
  code.evolve_model(tend)
  t2=time.time()

  print "time:",t2-t1

  print orbital_elements_from_binary(code.particles[0:2])
  
  print code.orbiters.position
  
  a,eps=elements(sun.mass,code.orbiters.x,code.orbiters.y,code.orbiters.z,
                     code.orbiters.vx,code.orbiters.vy,code.orbiters.vz,G=nbody_system.G)
  
  da=abs((a-a0)/a0)
  deps=abs(eps-eps0)/eps0

  print da,deps


  
if __name__=="__main__":
  for method in [1,0]:
    
    test_kepler_parabolic(tend=5000,method=method)
    print

#    test_kepler(N=1000000,tend=100.| units.yr,method=method)  

#  for method in [0,1]:
#    test_kepler(N=100000,tend=1000.| units.yr,method=method)  


#  for method in [0,1]:
#    test_kepler(N=100,tend=1e6| units.yr,method=method)  
  
