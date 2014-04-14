from amuse.units import units,constants,nbody_system

from amuse.community.mercury.interface import Mercury, MercuryWayWard
from amuse.community.huayno.interface import Huayno

from amuse.ext.orbital_elements import new_binary_from_orbital_elements,orbital_elements_from_binary

import numpy

from matplotlib import pyplot

mstar=1. | units.MSun
mdot=5.e-5 | units.MSun/units.yr

tend=15000.| units.yr

def starmass(t):
  return mstar - t*mdot

m1=0.001| units.MSun
m2=0.001| units.MSun

a1=200 | units.AU
a2=1000 | units.AU

e1=0.01
e2=0.5

be=new_binary_from_orbital_elements(mstar,m1,a1,eccentricity=e1, G=constants.G)
be[0].radius=1.| units.RSun
be[1].radius=6000. | units.km

mu1=constants.G*(mstar+m1)

tp=2*numpy.pi/mu1**(0.5)*a1**(1.5)
print "period:",tp.in_(units.day)
print "tend/period", tend/tp

be2=new_binary_from_orbital_elements(mstar+m1,m2,a2,eccentricity=0.5, G=constants.G)
be2[0].radius=1.| units.RSun
be2[1].radius=6000. | units.km

#be.add_particle(be2[1])

n=10000
dt=tend/n
dtplot=100*dt

code=Mercury()
code.parameters.timestep=dt/500
print code.parameters

code.particles.add_particles(be)


print dt/tp

t=code.model_time

pyplot.ion()
pyplot.show()

xa=numpy.arange(101)/100.
ya=10*7.66/(7.66-(7.66-1.438)*xa)
pyplot.plot(xa,ya)

while t<tend-dtplot/2:
  tnow=t
  while t< tnow+dtplot-dt/2:
    code.evolve_model(t+dt/2)
    code.particles[0].mass=starmass(t+dt)
    code.evolve_model(t+dt)
    t=code.model_time
    
  oe=orbital_elements_from_binary(code.particles[0:2], G=constants.G)
  print t/tend,oe[].in_(units.AU),oe[2]/(mstar/starmass(t)*a1),oe[3]
#  pyplot.plot(code.particles.x.value_in(units.AU),code.particles.y.value_in(units.AU),'r+')
#  pyplot.plot(no*t/tend,oe[2]/a,'r+')
#  pyplot.xlim(0,100)
#  pyplot.ylim(0.9,3.)
  pyplot.plot(t/tend,oe[3],'r+')
  pyplot.xlim(0,1.)
  pyplot.ylim(0,1.)
  pyplot.draw()
#  be=code.particles.copy()
#  dpos=be[1].position-be[0].position
#  dvel=be[1].velocity-be[0].velocity
#  print elements(be.mass.sum(),dpos[0],dpos[1],dpos[2],dvel[0],dvel[1],dvel[2])[1]

pyplot.savefig("test_2planet.png")

raw_input()

