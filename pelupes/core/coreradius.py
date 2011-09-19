import os
import sys
import numpy

from amuse.legacy.hop.interface import HopInterface as Hop

from matplotlib import pyplot

import logging
#logging.basicConfig(level=logging.DEBUG)

from amuse.ic.plummer import new_plummer_sphere
numpy.random.seed(1234)

def hopfromnb(nb,ids):
  m,r,x,y,z,vx,vy,vz,err=nb.get_state(ids)
  hop=Hop()
  ids2,err=hop.new_particle(x,y,z)    
  return hop,ids2

def plummer(x):
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

  tm=numpy.sum(mass)
  cmx=numpy.sum(mass*x)/tm
  cmy=numpy.sum(mass*y)/tm
  cmz=numpy.sum(mass*z)/tm

  return mass,radius,x,y,z,vx,vy,vz

def coreradius(mass,x,y,z):
  hop=Hop()
  ids,err=hop.new_particle(mass,x,y,z)
  hop.set_density_method(2)
  hop.set_nDens(7)
  hop.calculate_densities()
  dens,err=hop.get_density(ids)


  tdens=numpy.sum(dens)
  x_core=numpy.sum(dens*x)/tdens
  y_core=numpy.sum(dens*y)/tdens
  z_core=numpy.sum(dens*z)/tdens

#  r=((x-x_core)**2+(y-y_core)**2+(z-z_core)**2)**0.5
#  pyplot.figure(figsize=(8,6))
#  pyplot.loglog(r,dens,'r .')
#  pyplot.xlabel('r')
#  pyplot.ylabel('dens')
#  pyplot.show()
  
  rc=numpy.sqrt(
      numpy.sum(dens**2*((x-x_core)**2+(y-y_core)**2+(z-z_core)**2))/numpy.sum(dens**2))
  return x_core,y_core,z_core,rc
 

if __name__=="__main__":
  mass,radius,x,y,z,vx,vy,vz=plummer(100000)
  print coreradius(mass,x,y,z)
