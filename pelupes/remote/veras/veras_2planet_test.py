import numpy

from massloss import grav_with_se
from remote_se import forwarding_class_client
    
from amuse.units import units,nbody_system,constants
from amuse.datamodel import Particle

from amuse.ext.orbital_elements import new_binary_from_orbital_elements,orbital_elements_from_binary

from amuse.community.sse.interface import SSE

from matplotlib import pyplot

def veras_2planet_AGB_system(  a1=10 | units.AU):

  zams_mass=8. | units.MSun

  p=Particle(mass=zams_mass)
  
  se=SSE()
  p=se.particles.add_particle(p)
  
  se.evolve_model(41.9| units.Myr)
  
  m1=0.001| units.MSun
  m2=0.001| units.MSun
  
  a2=30 | units.AU
  
  e1=0.0
  e2=0.5
  
  be=new_binary_from_orbital_elements(p.mass,m1,a1,eccentricity=e1, G=constants.G)
  be[0].radius=p.radius
  be[1].radius=6000. | units.km
  
  mu1=constants.G*(be.total_mass())
    
  tp=2*numpy.pi/mu1**(0.5)*a1**(1.5)
  
  be2=new_binary_from_orbital_elements(be.total_mass(),m2,a2,eccentricity=0.5, G=constants.G)
  be2[1].radius=6000. | units.km
  
  be.add_particle(be2[1])
  
  be[0].is_star=True
  be[0].zams_mass=zams_mass
  
  return p.mass,p.age,be,tp
  
if __name__=="__main__":
  a1=10 | units.AU

  mstar,age,sys,tp=veras_2planet_AGB_system(a1)
    
  conv=nbody_system.nbody_to_si(1. | units.MSun,10.| units.AU)  
  
  timestep=tp*10
  dtplot=100*timestep
  dtend=(0.7 | units.Myr)
  tend=age+ dtend
  tnow=age

  print "period:",tp.in_(units.day)
  print "tend/period", dtend/tp
  print "tend/dtplot", dtend/dtplot
  
#  code=grav_with_se(converter=conv,timestep=timestep,begin_time=age)
  code=forwarding_class_client(grav_with_se,code_kwarg=dict(converter=conv,timestep=timestep,begin_time=age))

  print code.parameters["SSE"]
  print code.parameters["Huayno"]
    
  code.particles.add_particles(sys)
  
  tnow=code.model_time
  
  pyplot.ion()
  pyplot.show()
  
  xa=numpy.arange(101)/100.
  ya=10*7.66/(7.66-(7.66-1.438)*xa)
  pyplot.plot(xa,ya)
    
  p1=code.particles[[0,1]]
  p2=code.particles[[0,2]]  
  while tnow<tend-dtplot/2:
    code.evolve_model(tnow+dtplot)
    tnow=code.model_time
    
    print code.model_time.in_(units.Myr),
      
      
    t=tnow-age  
    oe=orbital_elements_from_binary(p1.get_intersecting_subset_in(code.particles), G=constants.G)
    print t/dtend,oe[2].in_(units.AU),(mstar/code.particles[0].mass)**-1,oe[3]

    pyplot.plot(t/dtend,oe[2].value_in(units.AU),'r+')
    oe=orbital_elements_from_binary(p2.get_intersecting_subset_in(code.particles), G=constants.G)
    pyplot.plot(t/dtend,oe[2].value_in(units.AU),'g+')
    pyplot.xlim(0,1.)
    pyplot.ylim(0,150.)
    pyplot.draw()


  
  pyplot.savefig("test_2planet.png")
  
  raw_input()

  
  
