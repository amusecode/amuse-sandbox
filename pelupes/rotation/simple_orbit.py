from amuse.couple.bridge import Bridge
from amuse.ext.bridge import bridge
from rotating_bridge import RotatingBridge
from amuse.datamodel import Particles,Particle

import numpy

from numpy import cos,sin

from amuse.units import nbody_system

import time

from matplotlib import pyplot

from amuse.ext.symplectic_composition import LEAPFROG,SPLIT_4TH_S_M6,SPLIT_8TH_SS_M21

class non_interacting_particles(object):
  def __init__(self,converter=None):
    self.particles=Particles()
    if converter is None:
      self.model_time=0| nbody_system.time
    else:
      self.model_time=converter.to_si(0 | nbody_system.time)  
  def evolve_model(self,tend):
    dt=tend-self.model_time
    self.particles.x+=dt*self.particles.vx
    self.particles.y+=dt*self.particles.vy
    self.particles.z+=dt*self.particles.vz
    self.model_time=tend
  @property
  def kinetic_energy(self):
    return 0.5*(self.particles.mass*(self.particles.vx**2+self.particles.vy**2+self.particles.vz**2)).sum()  
  @property
  def potential_energy(self):
    return 0. | (self.particles.mass.unit*self.particles.vx.unit**2)
    
class particle_potential(object):
  def __init__(self,M=1|nbody_system.mass, _G=nbody_system.G):
    self.M=M
    self._G=_G
  def get_gravity_at_point(self,eps,x,y,z):
    r2=x**2+y**2+z**2
    r=r2**0.5
    m=self.M  
    fr=self._G*m/(r2*r)
    ax=-fr*x
    ay=-fr*y
    az=-fr*z
    return ax,ay,az
  def get_potential_at_point(self,eps,x,y,z):
    r2=x**2+y**2+z**2
    r=r2**0.5
    m=self.M
    phi=-self._G*m/r
    return phi    
  def vcirc(self,r):  
    m=self.M  
    vc=(self._G*m/r)**0.5
    return vc

def inertial_to_rotating(t,omega,parts):
  x=parts.x
  y=parts.y
  vx=parts.vx
  vy=parts.vy
  rotating=parts.copy()
  rotating.x=x*cos(omega*t)+y*sin(omega*t)
  rotating.y=-x*sin(omega*t)+y*cos(omega*t)
  rotating.vx=(vx+y*omega)*cos(omega*t)+(vy-x*omega)*sin(omega*t)
  rotating.vy=-(vx+y*omega)*sin(omega*t)+(vy-x*omega)*cos(omega*t)
  return rotating
  
def rotating_to_inertial(t,omega,parts):     
  return inertial_to_rotating(t,-omega,parts)
   
def run(N,method=LEAPFROG):
   
  pot=particle_potential()   
  p=Particles(N)
  
  mu=pot.M*pot._G
  
  nstep=1000
  
  rmax=1 | nbody_system.length
  eps=0.9
  rmin=(1-eps)/(1+eps)*rmax
  a=(rmin+rmax)/2
  v0=mu**0.5*(2./rmax-1/a)**0.5
  T=2*numpy.pi*(a**3/mu)**0.5

  tnow=0 | nbody_system.time
  tend=T
  dt=T/nstep
  
  omega=(mu/a**3)**0.5
  
  p.mass=1 | nbody_system.mass
  p.x=rmax
  p.y=0*p.x
  p.z=0*p.x
  p.vy=v0
  p.vx=0*p.vy
  p.vz=0*p.vx
  p.radius=0|nbody_system.length
  
  rp=inertial_to_rotating(tnow,omega,p)
  
  ni=non_interacting_particles()
  p=ni.particles.add_particles(rp)
     
  sys=RotatingBridge(omega,method=method,use_threading=False)
  sys.add_system(ni, (pot,))
    
  
  sys.timestep=dt
  
  time=[tnow.number]
  x_r=[sys.particles.x.number]
  y_r=[sys.particles.y.number]

  inertial=rotating_to_inertial(tnow,omega,sys.particles)

  x_i=[inertial.x.number]
  y_i=[inertial.y.number]
  de=[0.]
  
  E0=sys.kinetic_energy+sys.potential_energy+sys.jacobi_potential_energy
  
  
  i=0
  while tnow<tend-dt/2:
    i+=1
    tnow+=dt
    sys.evolve_model(tnow)
    E=sys.kinetic_energy+sys.potential_energy+sys.jacobi_potential_energy
    
    time.append(tnow.number)
    
    x_r.append(sys.particles.x.number)
    y_r.append(sys.particles.y.number)

    inertial=rotating_to_inertial(tnow,omega,sys.particles)
    
    x_i.append(inertial.x.number)
    y_i.append(inertial.y.number)
    de.append(abs((E0-E)/E0))

  print i
  print max(de)
  
  pyplot.figure(figsize=(8,8))
  pyplot.subplot(211)
  pyplot.plot(x_i,y_i,'b-')
  pyplot.plot(x_r,y_r,'g-')
  pyplot.xlim(-2,2)
  pyplot.ylim(-2,2)
  pyplot.subplot(212)
  pyplot.semilogy(time,de)
  pyplot.show()


  
if __name__=="__main__":
#  import cProfile
#  cProfile.run("run()", "prof_full")
  run(1,method=LEAPFROG)
