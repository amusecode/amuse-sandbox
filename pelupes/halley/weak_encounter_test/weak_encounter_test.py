"""
this script tests the stabilizing effect of weak encounters

an ensemble of N Halleys (100 below) is integrated in the presence of a 
jupiter with mass x0, x0.2, x1.0 and x5.0 the real jupiter mass. The RMS 
spread in the resulting positions of the "Halleys" is recorded. 

In the resulting plot a linear growth of the RMS spread can be seen for 
the zero mass jupiter, hence no jupiter (as expected for gradual spread 
along the orbit due to initial conditions). The realistic jupiter mass 
shows a growth in the RMS consistent with the single simulation (shown 
in thepaper) i.e. linear or slow exponential growth in the beginning 
with an exponential growth at ~3000 yr

the sub jupiter mass simulation shows that a /smaller/ spread in the RMS 
compared with the no-jupiter sim. This shows that indeed the orbit of jupiter,
while chaotic in the presence of the current jupiter, would be stabalized if jupiter
where lower mass (transition mass must be somewhere between 0.2-1.0). This could
actually be detectable for either meteor swarms (dust connected to comets) or 
broken up comets/asteroids in certain orbits...

A High mass jupiter has a highle destabalizing effect on Halleys orbit (as expected)

suggestions: try the same with venus (shorter timescale), maybe 1 more jupiter mass (0.5?)

"""

import numpy

from amuse.units import units,nbody_system

from amuse.datamodel import Particles

from amuse.community.huayno.interface import Huayno

from matplotlib import pyplot

import cPickle

def soljup(mass_factor=1.):

  soljup=Particles(2)
  
  sol=soljup[0]
  jup=soljup[1]
  
  # 2457225.500000000 = A.D. 2015-Jul-22 00:00:00.0000 
  # data from horizons 
  # uncorrected for missing planets
  sol.mass=1.988544e30 | units.kg
  sol.x=0. | units.AU
  sol.y=0. | units.AU
  sol.z=0. | units.AU
  sol.vx=0. | units.AU/units.day
  sol.vy=0. | units.AU/units.day
  sol.vz=0. | units.AU/units.day
  sol.radius=6.955e5 | units.km
  
  jup.mass=1898.13e24 | units.kg
  jup.x=-4.674427698930987E+00 | units.AU
  jup.y=2.660151917962867E+00 | units.AU
  jup.z=9.354942380528984E-02 | units.AU
  jup.vx=-3.826060830850348E-03 | units.AU/units.day
  jup.vy=-6.207953616829373E-03 | units.AU/units.day
  jup.vz=1.113857181976940E-04 | units.AU/units.day
  jup.radius=69911. | units.km

  jup.mass*=mass_factor

  return soljup
  
def halleys(N=1,dp=0.,seed=123456):
  
  x=-2.047632191626166E+01 | units.AU
  y=2.540504375512831E+01 | units.AU
  z=-9.818400191369697E+00 | units.AU
  vx=-2.349849520378704E-05 | units.AU/units.day
  vy=8.785510030396358E-04 | units.AU/units.day
  vz=-1.524122463735410E-04 | units.AU/units.day
 
  p=Particles(N)
  p.mass=0. | units.MSun
  p.radius=0. | units.km
  p.x=x
  p.y=y
  p.z=z
  p.vx=vx
  p.vy=vy
  p.vz=vz
  
  r=(x**2+y**2+z**2)**0.5
  
  numpy.random.seed(seed)
  dx=r*numpy.random.normal(scale=dp,size=(N,3))
  
  p.position+=dx
  
  return p

def run(mass_factor=1.):


  conv=nbody_system.nbody_to_si(1 | units.MSun, 1. | units.AU)

  code=Huayno(conv)

  code.parameters.timestep_parameter=0.01
  print code.parameters

  ss=code.particles.add_particles(soljup(mass_factor=mass_factor))
  hs=code.particles.add_particles(halleys(N=100,dp=1.e-6))
  
  code.particles.move_to_center()  
  
  tend=10000. | units.yr
  dtplot= 5. | units.yr
  tnow=0. | units.yr
    
  t=[] | units.yr
  rms=[] | units.AU

  filename="t_rms_%4.2f"%mass_factor
  
  while tnow<tend:
    code.evolve_model(tnow+dtplot)
    tnow=code.model_time
    print tnow/tend
    #~ pyplot.clf()
    #~ pyplot.subplot(121)
    #~ pyplot.plot(ss.x.value_in(units.AU),ss.y.value_in(units.AU),'r+')
    #~ pyplot.plot(hs.x.value_in(units.AU),hs.y.value_in(units.AU),'b.')
    #~ pyplot.xlim(-30,30)
    #~ pyplot.ylim(-30,30)
    rms.append(numpy.std(hs.position,axis=0).length())
    t.append(tnow)
    #~ pyplot.subplot(122)
    #~ pyplot.semilogy(t/tend,rms.value_in(units.AU))
    #~ pyplot.xlim(0,1)
    #~ pyplot.ylim(1.e-6,10.)
    #~ pyplot.draw()
#~ 
  f=open(filename,'w')
  cPickle.dump([t,rms],f)
  f.close()
  
if __name__=="__main__":
  pyplot.ion()
  f=pyplot.figure(figsize=(12,6))
  pyplot.show()

  run(0.5)
  #~ run(0.2)
  #~ run(1.)
  #~ run(5.)
