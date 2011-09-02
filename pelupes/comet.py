import numpy
import cPickle

import bridge
import tidal_field


from amuse.units import nbody_system as NU
from amuse.community.twobody import twobody as interface

from matplotlib import pyplot

from amuse.units import *
convert1 = NU.nbody_to_si(10**9 | MSun, 1000 | parsec)
convert2 = NU.nbody_to_si(1 | MSun, 100000 | AU)

convert_1_2=lambda x : convert2.to_nbody(convert1.to_si(x)) 

def total_energy_perturbed(nb, perturber):
  state,err=nb.get_state(0)
  time,err=nb.get_time()
  Ep,err=nb.get_potential_energy()
  Ek,err=nb.get_kinetic_energy()
  result,err=perturber.get_potential_at_point(0.,state['x'],state['y'],state['z'])  
  Es=result
  Etot=Ep+Es+Ek
  print "time = %.5f, Ep = %.6f, Ek = %.6f, Ep = %.6f, Etot = %.6f" % \
      (time,Ep,Ek,Es,Etot)

def tidal_field_from_file(fnaam):
  f=open(fnaam,'r')
  time=convert_1_2( (NU.time).new_quantity(cPickle.load(f)) )
  tides=convert_1_2( (NU.time**-2).new_quantity(cPickle.load(f)) )
  f.close() 
  tf=tidal_field.time_dependent_tidal_field(time.number,tides.number)
  return tf

def make_plot(x,y):
  pyplot.figure(figsize=(8,8))
  pyplot.plot(x,y)
  pyplot.plot([0.],[0.],'x')
  pyplot.axis([-1.2,1.2,-0.5,1.9])
  pyplot.xlabel('( x100k AU )')
  pyplot.savefig('test.png')  

def tide_test():
  nb = interface.twobody()
  nb.setup_module()
  
  mu=1.
  rinit=.2
  vinit=2.95

  energy=0.5*vinit**2-mu/rinit
  a=-mu/(2*energy)
  ecc=1-rinit/a
  rmax=a*(1+ecc)
  print a,ecc,rmax
  tperiod=2*numpy.pi/numpy.sqrt(mu)*a**(3./2)

  tidal=tidal_field_from_file('tideseq-1013801.dat')
  
  nb.new_particle( 1.,0.,rinit,0.,0.,0.,vinit,0.)

  dt=.001*tperiod
  t=0.

  total_energy_perturbed(nb, tidal)
  
  
  x=[]
  y=[]
  while t < tperiod*10:
    t=t+dt 
    bridge.evolve_w_kick(t,nb,tidal,ids=[0])
    state,err=nb.get_state(0)
    x.append( state['y'])
    y.append( -state['x'])
    
  make_plot(x,y)  
    
  print convert2.to_si(t|NU.time).in_(Myr)

  total_energy_perturbed(nb, tidal)


if __name__=="__main__":
  tide_test()
