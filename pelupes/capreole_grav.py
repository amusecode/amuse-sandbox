import os
import sys

from amuse.community.capreole.interface import CapreoleInterface
from amuse.community.capreole.interface import GLCapreoleInterface

import numpy

from matplotlib import pyplot

gamma=5./3
ggrav=980.

def constantdens_constantgrav(N,L):
  instance=GLCapreoleInterface()
  instance.initialize_code()
  instance.setup_mesh(N/10,N/10,N,L/10.,L/10.,L)
  instance.commit_parameters()
  x,y,z=numpy.indices( (N/10,N/10,N) )
  x=x.flatten()+1
  y=y.flatten()+1
  z=z.flatten()+1
  rho=0.001*numpy.ones_like(x)
  rhvx=0.*numpy.ones_like(x)
  rhvy=0.*numpy.ones_like(x)
  rhvz=0.*numpy.ones_like(x)
  en=1.e6*numpy.ones_like(x)
  instance.set_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
  
  gx=0.*numpy.ones_like(x)
  gy=0.*numpy.ones_like(x)
  gz=-980.*numpy.ones_like(x)
  err=instance.set_gravity_field(x,y,z,gx,gy,gz)

  instance.initialize_grid()
  instance.viewer()
  return instance

def eqdens_constantgrav(N,L):
  instance=GLCapreoleInterface(debugger='xterm')
  instance.initialize_code()
  instance.setup_mesh(N/10,N/10,N,L/10.,L/10.,L)
  instance.commit_parameters()
  ix,iy,iz=numpy.indices( (N/10,N/10,N) )
  ix=ix.flatten()+1
  iy=iy.flatten()+1
  iz=iz.flatten()+1
  x,y,z,err=instance.get_position_of_index(ix,iy,iz)
  rho0=0.001
  pres0=1.e6
  c2=gamma*pres0/rho0
  zscl=c2/ggrav/gamma
  rho=rho0*numpy.exp(-z/zscl)
  rhvx=0.*numpy.ones_like(rho)
  rhvy=0.*numpy.ones_like(rho)
  rhvz=0.*numpy.ones_like(rho)
  en=pres0*numpy.exp(-z/zscl)/(gamma-1.)
  
  instance.set_grid_state(ix,iy,iz,rho,rhvx,rhvy,rhvz,en)

  ix,iy,iz=numpy.indices( (N/10,N/10,N) )
  ix=ix.flatten()+1
  iy=iy.flatten()+1
  iz=iz.flatten()+1
  gx=0.*numpy.ones_like(ix)
  gy=0.*numpy.ones_like(ix)
  gz=-ggrav*numpy.ones_like(ix)
  err=instance.set_gravity_field(ix,iy,iz,gx,gy,gz) 
  rs="ref_shift"
  rf="reflective"
  err=instance.set_boundary(rf,rf,rf,rf,rf,rf)
  print err
  instance.initialize_grid()
  if hasattr(instance,'viewer'):
    instance.viewer()
  return instance


def plot(x,rho,pres,vel,U):
  pyplot.subplot( 221)
  pyplot.plot(x,rho)
  pyplot.subplot( 222)
  pyplot.plot(x,pres)
  pyplot.subplot( 223)
  pyplot.plot(x,vel)
  pyplot.subplot( 224)
  pyplot.plot(x,u)
  pyplot.savefig("plot.png")

if __name__=="__main__":
  N=100
  instance=eqdens_constantgrav(N,1.e7)
  instance.evolve_model(1000.)
  ix,iy,iz=numpy.indices( (1,1,100) )
  ix=ix.flatten()+N/20  
  iy=iy.flatten()+N/20  
  iz=iz.flatten()+1
  x,y,z,err=instance.get_position_of_index(ix,iy,iz)
  rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(ix,iy,iz)
  vel=rhovz/rho
  pres=(gamma-1)*(en-0.5*rho*vel**2)
  u=pres/(gamma-1)/rho
  plot(z/1.e5,rho,pres,vel,u)
  instance.stop()
