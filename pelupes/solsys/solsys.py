import numpy
import cPickle
from time import time
from amuse.units import units
from amuse.units import nbody_system

from amuse.ext.solarsystem import new_solar_system,new_solar_system_for_mercury

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

from amuse.community.adaptb.interface import Adaptb
from amuse.community.bhtree.interface import BHTree
from amuse.community.huayno.interface import Huayno
from amuse.community.ph4.interface import ph4
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.mercury.interface import MercuryWayWard,Mercury

def solsys(interface, tend=1000. | units.yr, dt=1 | units.day, parameters=[]):

  sss = new_solar_system()

  convert=nbody_system.nbody_to_si( 1.| units.MSun, 5.| units.AU)
  grav=interface(convert,channel_type="mpi")
  grav.initialize_code()

  for name,value in parameters:
    setattr(grav.parameters, name, value)
  
  grav.particles.add_particles(sss)
  grav.commit_particles()

  try:
    E0=grav.kinetic_energy+grav.potential_energy
  except:
    E0=grav.particles.kinetic_energy()+grav.particles.potential_energy()

  E=[1.e-20]
  time=[0.]
  
  x=grav.particles.x.value_in(units.AU)
  xx=[x]
  y=grav.particles.y.value_in(units.AU)
  yy=[y]
  
  t=0 | units.yr
  dt=convert.to_si(dt).in_(units.yr)
  
  i=0
  while t<tend-dt/2:
    i+=1
    t+=dt
    grav.evolve_model(t)
    if i%100==0:
      print (t/tend)
    
    try:
      e=grav.kinetic_energy+grav.potential_energy
    except:
      e=grav.particles.kinetic_energy()+grav.particles.potential_energy()
    E.append( abs((e-E0)/E0) )
    time.append( t.value_in(units.yr) )
    x=grav.particles.x.value_in(units.AU)
    xx.append(x)
    y=grav.particles.y.value_in(units.AU)
    yy.append(y)
  
  grav.stop()  
    
  f=pyplot.figure( figsize=(8,6))  
  time=numpy.array(time)
  E=numpy.array(E)
    
  xx=numpy.array(xx)
  yy=numpy.array(yy)
  
  return time,E,xx,yy

if __name__=="__main__":
  import cProfile
  time,E,xx,yy=solsys(Hermite, tend=10. | units.yr, dt=1. | units.day,parameters=[("timestep_parameter",0.02)])
#  cProfile.run("solsys(Hermite, tend=10. | units.yr, dt=1. | units.day,parameters=[('timestep_parameter',0.02)])","prof")
