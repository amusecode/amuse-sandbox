import numpy

from amuse.community.twobody import twobody
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse import datamodel
def binary(m1=1.|units.MSun, m2=1|units.MSun, r1=1.| units.RSun, r2=1.|units.RSun, \
             period=1.| units.yr,ecc=0):
  convert_nbody = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.AU)
    
  mu=constants.G*(m1+m2)
# semi major axis
  a=(period**2/(4*numpy.pi**2)*mu)**(1./3.)  
  print "semi major axes:", a.value_in(units.AU)
  
# aphelion
  ra=a*(1+ecc)
  
# angular momentum
  h=(mu*a*(1-ecc**2))**0.5  
  
# vel at aphelion
  va=h/ra  
  
  x1=m2/(m1+m2)*ra
  y1=z1=0.*x1
  vy1=-m2/(m1+m2)*va
  vx1=vz1=0.*vy1

  x2=-m1/(m1+m2)*ra
  y2=z2=0.*x2
  vy2=m1/(m1+m2)*va
  vx2=vz2=0.*vy2
  
  parts=datamodel.Particles(2)

  parts[0].mass=m1
  parts[0].radius=r1
  parts[0].x=x1
  parts[0].y=y1
  parts[0].z=z1
  parts[0].vx=vx1
  parts[0].vy=vy1
  parts[0].vz=vz1

  parts[1].mass=m2
  parts[1].radius=r2
  parts[1].x=x2
  parts[1].y=y2
  parts[1].z=z2
  parts[1].vx=vx2
  parts[1].vy=vy2
  parts[1].vz=vz2
 
  gravity = twobody.TwoBody(convert_nbody)
 
  gravity.particles.add_particles(parts)
  
  print gravity.particles
  gravity.evolve_model(period/2.)
  print gravity.particles
  gravity.evolve_model(period)
  print gravity.particles

if __name__=="__main__":
  binary()
