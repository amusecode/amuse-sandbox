import os
import sys
import numpy
import cPickle


from amuse.ic.plummer import MakePlummerModel
from amuse.ic.plummer import new_plummer_model
from amuse.community.phiGRAPE.interface import PhiGRAPEInterface as phi
from amuse.community.hop.interface import HopInterface as Hop

#from matplotlib import pyplot

import logging
#logging.basicConfig(level=logging.DEBUG)


numpy.random.seed(123456)

def get_lagrangian_radii(grav,ids):
  m,r,x,y,z,vx,vy,vz,err=grav.get_state(ids)
  tm=numpy.sum(m)
  cmx=numpy.sum(m*x)/tm
  cmy=numpy.sum(m*y)/tm
  cmz=numpy.sum(m*z)/tm
  print 'cmx1',cmx,cmy,cmz
  r2=(x-cmx)**2+(y-cmy)**2+(z-cmz)**2
  rsorted=numpy.sqrt(numpy.sort(r2))
  r01=rsorted[long(len(x)*0.01)]
  r02=rsorted[long(len(x)*0.02)]
  r05=rsorted[long(len(x)*0.05)]
  r10=rsorted[long(len(x)*0.1)]
  r20=rsorted[long(len(x)*0.2)]
  r50=rsorted[long(len(x)*0.5)]
  r75=rsorted[long(len(x)*0.75)]
  print "lr:",r01,r02,r05,r10,r20,r50,r75
  r2=(x-cmx)**2+(y-cmy)**2+(z-cmz)**2
  r2sorted=numpy.sort(r2)
  rm=r2sorted[long(len(x)*0.90)]
  a=r2 < rm
  mm=m.compress(a)
  xx=x.compress(a)
  yy=y.compress(a)
  zz=z.compress(a)
  tmm=numpy.sum(mm)
  cmx=numpy.sum(mm*xx)/tmm
  cmy=numpy.sum(mm*yy)/tmm
  cmz=numpy.sum(mm*zz)/tmm
  print 'cmx2',cmx,cmy,cmz
  r2=(x-cmx)**2+(y-cmy)**2+(z-cmz)**2
  rsorted=numpy.sqrt(numpy.sort(r2))
  r01=rsorted[long(len(x)*0.01)]
  r02=rsorted[long(len(x)*0.02)]
  r05=rsorted[long(len(x)*0.05)]
  r10=rsorted[long(len(x)*0.1)]
  r20=rsorted[long(len(x)*0.2)]
  r50=rsorted[long(len(x)*0.5)]
  r75=rsorted[long(len(x)*0.75)]
  print "lr2:",r01,r02,r05,r10,r20,r50,r75
  return [r01,r02,r05,r10,r20,r50,r75]

def hopfromnb(nb,ids):
  m,r,x,y,z,vx,vy,vz,err=nb.get_state(ids)
  hop=Hop()
  ids2,err=hop.new_particle(x,y,z)    
  return hop,ids2

def centre(nb,ids):
  hop,ids2=hopfromnb(nb,ids)
  hop.calculate_densities()
  dens,err=hop.get_density(ids2)
  mi=dens.argmax()
  x,y,z,err=hop.get_position(ids2[mi])
  return x,y,z

def plummer(x,interface):
  plummer=MakePlummerModel(x)
  mass,pos,vel=plummer.new_model()

  mass=mass[0:,0]
  x=pos[0:,0]
  y=pos[0:,1]
  z=pos[0:,2]

  vx=vel[0:,0]
  vy=vel[0:,1]
  vz=vel[0:,2]
  radius=mass*0.
  nb = interface(redirection="none")#,debugger="gdb")
  nb.initialize_code()

  ids,error = nb.new_particle(mass,radius,x,y,z,vx,vy,vz)
  if filter(lambda x: x != 0, error) != []: raise Exception

  return nb,ids

if __name__=="__main__":
  import time
  nb,ids=plummer(1024,phi)
#  nb.set_inttype_parameter(4)
#  nb.set_timestep_parameter(0.01)
  nb.set_eps2(1.e-6)
  nb.commit_particles()
  ek,err=nb.get_kinetic_energy()
  ep,err=nb.get_potential_energy()
  e1=ek+ep
  print ek,ep,ek+ep
  t1=time.time()

  print centre(nb,ids)
