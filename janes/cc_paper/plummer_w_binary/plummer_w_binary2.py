# test of huayno with fractal ic

"""
Energy error and runtime as a function of fractal dimension for 
different integrators. Plotted are the energy error (top panel) and
runtime (bottom panel) of runs with 1024 particles as a 
function of the fractal dimension of the particle distribution for the 
EXTRAPOLATE (dashed line), SF-split (dash-dotted line) and the CC 
(dotted) and CC_KEPLER (drawn) split integrators. Note that all 
integrators have similar error behaviour, while the run time increases 
for decreasing fractal dimension (so for more structured particle 
distributions) for the EXTRAPOLATE and SF split integrators, while the 
runtime is flat and even decreases slightly for lower fractal dimension 
for the CC integrators. 

"""

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
from amuse.io import write_set_to_file    
    
from amuse.community.fi.interface import Fi
from amuse.community.huayno.interface import Huayno
from amuse.community.ph4.interface import ph4
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE

from amuse.ic.plummer import new_plummer_model
from amuse.ext.orbital_elements import new_binary_from_orbital_elements

from amuse.datamodel import Particles

from time import time,sleep

from numpy import random

def plummer_w_binaries(N, f_bin=1.,logamin=-4,logamax=-0.5,seed=123454321):
  random.seed(seed)
  plummer=new_plummer_model(N)
  semi=10**random.uniform(logamin,logamax,N) | nbody_system.length
  inc=numpy.arccos(random.uniform(-1.,1.,N))/numpy.pi*180
  longi=random.uniform(0,2*numpy.pi,N)/numpy.pi*180
  binaries=Particles()
  tobin=plummer[:int(f_bin*N)]
  nobin=plummer-tobin
  if len(nobin)>0:
    binaries.add_particles(nobin)  
  for i,p in enumerate( tobin ):
    mass1=p.mass/2
    mass2=p.mass/2
    binary=new_binary_from_orbital_elements(mass1,mass2,semi[i],
              inclination=inc[i],longitude_of_the_ascending_node=longi[i])
    binary.position+=p.position
    binary.velocity+=p.velocity
    binaries.add_particles(binary)
  return binaries  


def timerun(N=1000,f_bin=1.,tend=0.25| nbody_system.time,inttype=Huayno.inttypes.HOLD_DKD,seed=12345678,logamin=-4):

  stars = plummer_w_binaries(N=N,f_bin=f_bin,seed=seed,logamin=logamin)
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
  return t2-t1,abs((E0-E)/E0),p0-p

def generate_data(Nsample=10,N=128,f_bin=1., label=""):
  data=dict()
  data['N']=N
  logamins=[-0.5,-1.,-1.5,-2.,-2.5,-3.]#,-3.5,-4.]
  data['logamins']=logamins
  inttypes=[Huayno.inttypes.HOLD_DKD,Huayno.inttypes.CC,Huayno.inttypes.CC_KEPLER,Huayno.inttypes.EXTRAPOLATE]
#  inttypes=[Huayno.inttypes.HOLD_DKD]
  data['inttypes']=inttypes
  for inttype in inttypes:
    data[inttype]=dict()
    for logamin in logamins:
      print logamin
      data[inttype][logamin]=[]
      for i in range(Nsample):
        seed=123456+i
        t,de,dp=timerun(N=N,f_bin=f_bin,inttype=inttype,seed=seed,logamin=logamin)
        print t,de,dp.length()
        data[inttype][logamin].append( (t,de,dp) )
    print     
  f=open("data"+label,"wb")
  cPickle.dump(data,f)
  f.close()      

def doplot(label=""):
  f=open("data"+label,"rb")
  data=cPickle.load(f)
  f.close()
  inttypes=data['inttypes']
  for inttype in inttypes:
    N=data['N']
    logamins=data['logamins']
    time=[]
    time_disp=[]
    eerr=[]
    eerr_disp=[]
    minerr=[]
    maxerr=[]
    for logamin in logamins:
      t=numpy.array(map(lambda x: x[0],data[inttype][logamin]))
      de=numpy.array(map(lambda x: abs(x[1]),data[inttype][logamin]))
      time.append(numpy.average(t))
      time_disp.append(numpy.std(t))
      eerr.append(numpy.average(de))
      minerr.append(numpy.min(de))
      maxerr.append(numpy.max(de))
      eerr_disp.append(numpy.std(de))
    data[inttype]["logamins"]=numpy.array(logamins)  
    data[inttype]["time"]=numpy.array(time)  
    data[inttype]["timedisp"]=numpy.array(time_disp)
    data[inttype]["eerr"]=numpy.array(eerr)
    data[inttype]["eerrdisp"]=numpy.array(eerr_disp)
    data[inttype]["minerr"]=numpy.array(minerr)
    data[inttype]["maxerr"]=numpy.array(maxerr)
  
  linestyles=["-.","-",":","--"]
  colors="rgbcyk"
    
  f=pyplot.figure(figsize=(8,10))
  ax=f.add_subplot(211)
  for i,inttype in enumerate(inttypes):
    logamins=data[inttype]['logamins']
    eerr=data[inttype]['eerr']
    minerr=data[inttype]['minerr']
    maxerr=data[inttype]['maxerr']
    ax.loglog(10**logamins,eerr,'r'+linestyles[i])
#    ax.semilogy(fdims,minerr,colors[i]+":",lw=0.5)
#    ax.semilogy(fdims,maxerr,colors[i]+":",lw=0.5)
    ax.set_ylabel("|E-E0|/E0")
  ax=f.add_subplot(212)
  for i,inttype in enumerate(inttypes):
    logamins=data[inttype]['logamins']
    eerr=data[inttype]['time']
    ax.loglog(10**logamins,eerr,'r'+linestyles[i])
    ax.set_xlabel("log minimum a")
    ax.set_ylabel("wallclock time (s)")
    
  pyplot.savefig("plummer_w_binary"+label+'.eps')
    
if __name__=="__main__":
  generate_data(Nsample=3,N=64,label="128")
#  generate_data(Nsample=10,logamin=-3,label="-3")
  doplot(label="128")

# plot as a function of fd: dE , time_CC/ time_[HOLD_DKD,EXTRAPOLATE]
# (mean, distribution of 10 random seed?)

