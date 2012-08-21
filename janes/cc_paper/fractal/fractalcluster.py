# test of huayno with fractal ic

import os
import sys
import numpy
import cPickle

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    
from amuse.units import nbody_system
from amuse.units import units
from amuse.io import write_set_to_file    
    
from amuse.community.fi.interface import Fi
from amuse.community.huayno.interface import Huayno
from amuse.community.ph4.interface import ph4
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE

from amuse.ic.fractalcluster import new_fractal_cluster_model

from time import time,sleep


def timefractal(N=1000,fd=1.6,tend=0.25| nbody_system.time,inttype=Huayno.inttypes.HOLD_DKD,seed=12345678):

  stars = new_fractal_cluster_model(N=N,fractal_dimension=fd,random_seed=seed,
                                    virial_ratio=0.5,do_scale=True,verbose=False)
  code=Huayno()
  code.parameters.inttype_parameter=inttype
  code.parameters.timestep_parameter=0.01
  code.particles.add_particles(stars)

  E0=code.kinetic_energy+code.potential_energy
  p0=code.particles.total_momentum()
  t1=time()
  code.evolve_model(tend)
  t2=time()
  E=code.kinetic_energy+code.potential_energy
  p=code.particles.total_momentum()
  code.stop()
  return t2-t1,abs(E0-E)/E0,p0-p

def generate_data():
  N=1024
  data=dict()
  for inttype in [Huayno.inttypes.HOLD_DKD,Huayno.inttypes.CC,Huayno.inttypes.CC_KEPLER,Huayno.inttypes.EXTRAPOLATE]:
    data[inttype]=dict()
    for fd in [1.6,2.0,2.3,2.6,2.8,3.0]:
      data[inttype][fd]=[]
      for i in range(25):
        seed=123456+i
        t,de,dp=timefractal(N=N,fd=fd,inttype=inttype,seed=seed)
        print t,de,dp.length()
        data[inttype][fd].append( (t,de,dp) )
    print     
  f=open("data","wb")
  cPickle.dump(data,f)
  f.close()      

def doplot():
  f=open("data","rb")
  data=cPickle.load(f)
  f.close()
  inttypes=[Huayno.inttypes.HOLD_DKD,Huayno.inttypes.CC_KEPLER,Huayno.inttypes.CC,Huayno.inttypes.EXTRAPOLATE]
  for inttype in inttypes:
    fdims=[1.6,2.0,2.3,2.6,2.8,3.0]
    time=[]
    time_disp=[]
    eerr=[]
    eerr_disp=[]
    for fd in fdims:
      t=numpy.array(map(lambda x: x[0],data[inttype][fd]))
      de=numpy.array(map(lambda x: abs(x[1]),data[inttype][fd]))
      time.append(numpy.average(t))
      time_disp.append(numpy.disp(t))
      eerr.append(numpy.average(de))
      eerr_disp.append(numpy.disp(de))
    data[inttype]["fdims"]=fdims  
    data[inttype]["time"]=time  
    data[inttype]["timedisp"]=time_disp
    data[inttype]["eerr"]=eerr
    data[inttype]["eerrdisp"]=eerr_disp  
  
  linestyles=["-.","-",":","--"]
    
  f=pyplot.figure(figsize=(8,10))
  ax=f.add_subplot(211)
  for i,inttype in enumerate(inttypes):
    fdims=data[inttype]['fdims']
    eerr=data[inttype]['eerr']
    ax.semilogy(fdims,eerr,'r'+linestyles[i])
    ax.set_ylabel("|E-E0|/E0")
  ax=f.add_subplot(212)
  for i,inttype in enumerate(inttypes):
    fdims=data[inttype]['fdims']
    eerr=data[inttype]['time']
    ax.semilogy(fdims,eerr,'r'+linestyles[i])
    ax.set_xlabel("fractal dimension")
    ax.set_ylabel("wallclock time (s)")
    
  pyplot.savefig('fractal_dimension_dE_time.eps')
    
if __name__=="__main__":
  generate_data()
  doplot()

# plot as a function of fd: dE , time_CC/ time_[HOLD_DKD,EXTRAPOLATE]
# (mean, distribution of 10 random seed?)

