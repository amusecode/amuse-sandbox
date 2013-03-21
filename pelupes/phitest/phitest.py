import os
import sys
import numpy
import cPickle
import time

from matplotlib import pyplot

from amuse.community.huayno.interface import Huayno
from amuse.community.ph4.interface import ph4
from amuse.community.hermite0.interface import Hermite

from amuse.units import nbody_system

from amuse.ic.plummer import new_plummer_model


import logging
#logging.basicConfig(level=logging.DEBUG)

numpy.random.seed(456789)

def diagnostics(data,particles,modeltime,wc,ek,ep):
    core=particles.cluster_core(reuse_hop=True)
    
    data.setdefault('mf',[0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 0.75, 0.9, 1])

    lr=particles.LagrangianRadii(cm=core.position,reuse_hop=True,mf=data['mf'])
    
    cm=particles.center_of_mass()
    cmv=particles.center_of_mass_velocity()
    totalp=particles.total_momentum()
    totall=particles.total_angular_momentum()
    
    data.setdefault('time',[]).append(modeltime)
    data.setdefault('kinetic_energy',[]).append(ek) 
    data.setdefault('potential_energy',[]).append(ep) 
    data.setdefault('core',[]).append(core.as_set())
    data.setdefault('lagrangian_radii',[]).append(lr[0])
    data.setdefault('center_of_mass',[]).append(cm)
    data.setdefault('center_of_mass_velocity',[]).append(cmv)
    data.setdefault('total_momentum',[]).append(totalp)
    data.setdefault('total_angular_momentum',[]).append(totall)
    data.setdefault('wallclock',[]).append(wc)

def run_to_core_collaps(code,N=1024,
                        tend=400 | nbody_system.time,
                        dt=0.25 | nbody_system.time,
                        startup=dict(),
                        parameters=dict()):
  
  plum=new_plummer_model(N)

  nb=code(**startup)  

  print nb.parameters

  for key,value in parameters.items():
    setattr(nb.parameters,key,value)

  nb.particles.add_particles(plum)

  ek=nb.kinetic_energy
  ep=nb.potential_energy

  print ek,ep

  e2=ek+ep

  deltaE=0.*e2
  Eshift=0.*e2
  E0=e2

  eta=0.
  
  t=0. | nbody_system.time
  i=0

  t2=time.time()
  t1=t2

  data=dict()
#  diagnostics(data,nb.particles.copy(),nb.model_time,t2-t1,ek,ep)
  while (t<tend-dt/2):
    i=i+1
    t=t+dt
    nb.evolve_model(t)

    ek=nb.kinetic_energy
    ep=nb.potential_energy
    
    e1=e2
    e2=ek+ep
    deltaE+=e2-e1
    print t,':',(e2-e1)/E0,deltaE/E0,(e2-Eshift-E0)/E0
    tp=t2
    t2=time.time()
    if eta==0.:
      eta=(t2-tp)/dt*(tend-t)
    eta=0.75*eta+0.25*((t2-tp)/dt*(tend-t))
    print 'ETA:', eta/3600.
    
#    diagnostics(data,nb.particles.copy(),nb.model_time,t2-t1,ek,ep)

  t2=time.time()

  f=open('data','wb')
  cPickle.dump(data,f)
  f.close()

  time.sleep(0.1)
  ek=nb.kinetic_energy
  ep=nb.potential_energy
  e2=ek+ep
  nb.stop()
  time.sleep(0.1)

  print 'total time:', t2-t1

#  execfile('./plot.py')


if __name__=="__main__":
  
#  print dir(Huayno.inttypes)
#  raise
  
#  run_to_core_collaps(Hermite,N=128,tend=30.| nbody_system.time,parameters=dict())
  run_to_core_collaps(Huayno,N=128,tend=30.| nbody_system.time,
     parameters=dict(timestep_parameter=0.01),startup=dict(channel_type="sockets"))
#  run_to_core_collaps(ph4,N=128,tend=30.| nbody_system.time,
#      parameters=dict(timestep_parameter=0.1))
