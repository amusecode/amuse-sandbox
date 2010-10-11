import os
import sys

from amuse.legacy.capreole.interface import CapreoleInterface
from amuse.legacy.capreole.interface import GLCapreoleInterface

import numpy
from matplotlib import pyplot
import matplotlib.cm as cm

def cloud_shock_problem(N=32,M=4,xi=10.,L=10.,rc=1.,gamma=5./3):

  cs=numpy.sqrt((gamma-1))
  cs_out=numpy.sqrt((gamma-1)*xi)
  vs=cs_out*2.7

  dl=L/N
  dr=numpy.sqrt(3.)*dl

  ix,iy,iz=numpy.indices( (N,4*N,N) )

  ix=ix.reshape((4*N*N*N,))+1
  iy=iy.reshape((4*N*N*N,))+1
  iz=iz.reshape((4*N*N*N,))+1

  x=L*(ix/(1.*N)-0.5-1/2./N)
  y=L*(iy/(1.*N)-0.5-1/2./N)
  z=L*(iz/(1.*N)-0.5-1/2./N)

  r=(x**2+y**2+z**2)**0.5
  rho=numpy.zeros_like(r)

  selection=r>(rc+dr)
  selection=numpy.compress(selection,numpy.arange(len(r))) 
  rho[selection]=1./xi
  selection=r<(rc-dr)
  selection=numpy.compress(selection,numpy.arange(len(r))) 
  rho[selection]=1.

  ix,iy,iz=numpy.indices( (M,M,M) )
  ix=ix.reshape((M*M*M,))+1
  iy=iy.reshape((M*M*M,))+1
  iz=iz.reshape((M*M*M,))+1

  dx=dl*(ix/(1.*M)-0.5-1/2./M)
  dy=dl*(iy/(1.*M)-0.5-1/2./M)
  dz=dl*(iz/(1.*M)-0.5-1/2./M)

  drho=numpy.zeros_like(dx)

  selection= numpy.logical_and( r>=(rc-dr), r <= (rc+dr))
  selection=numpy.compress(selection,numpy.arange(len(r))) 


  for i in selection:
    xx=x[i]  
    yy=y[i]
    zz=z[i]
    r=((xx+dx)**2+(yy+dy)**2+(zz+dz)**2)**0.5
    dselection=r>rc
    dselection=numpy.compress(dselection,numpy.arange(len(r))) 
    drho[dselection]=1./xi
    dselection=r<=rc
    dselection=numpy.compress(dselection,numpy.arange(len(r))) 
    drho[dselection]=1.
    rho[i]=drho.mean()

  rho=rho.reshape((N,4*N,N))

  rhvx=numpy.zeros_like(rho)
  rhvy=numpy.zeros_like(rho)
  rhvz=numpy.zeros_like(rho)
  rhvy[0:N,0:3,0:N]=1./xi*vs
  en=1./rho
  en[0:N,0:3,0:N]=en[0:N,0:3,0:N]+1./xi*vs**2

  ix,iy,iz=numpy.indices( (N,4*N,N) )

  ix=ix+1
  iy=iy+1
  iz=iz+1

  return ix,iy,iz,rho,rhvx,rhvy,rhvz,en

def cloudshock(N=32,M=4,xi=10.,L=10.,rc=1.,gamma=5./3):
  ix,iy,iz,rho,rhvx,rhvy,rhvz,en=cloud_shock_problem(N=N,M=M,xi=xi,L=L,rc=rc,gamma=gamma)
  ix=ix.reshape( (4*N*N*N,) )
  iy=iy.reshape( (4*N*N*N,) )
  iz=iz.reshape( (4*N*N*N,) )
  rho=rho.reshape( (4*N*N*N,) )
  rhvx=rhvx.reshape( (4*N*N*N,) )
  rhvy=rhvy.reshape( (4*N*N*N,) )
  rhvz=rhvz.reshape( (4*N*N*N,) )
  en=en.reshape( (4*N*N*N,) )
  instance=GLCapreoleInterface()
  instance.initialize_code()
  instance.set_boundary("reflective","reflective","outflow","outflow",
                          "reflective","reflective")
  instance.setup_mesh(N,4*N,N,L,4*L,L)
  instance.commit_parameters()
  instance.set_grid_state(ix,iy,iz,rho,rhvx,rhvy,rhvz,en)
  instance.initialize_grid()
  if hasattr(instance,'viewer'):
    instance.viewer()
  return ix,iy,iz,instance

if __name__=="__main__":
  N=32
  ix,iy,iz,instance=cloudshock(N=N)
  
  rc=1.
  gamma=5./3
  xi=10.
  cs=numpy.sqrt((gamma-1))
  cs_out=numpy.sqrt((gamma-1)*xi)
  vs=cs_out*2.7
  tau=1.6*2*rc*xi**0.5/vs
  
  ix=ix.reshape( (N,4*N,N) )
  iy=iy.reshape( (N,4*N,N) )
  iz=iz.reshape( (N,4*N,N) )
  
  ix_slice=ix[0:N,0:4*N,N/2].reshape( (N*4*N,) )
  iy_slice=iy[0:N,0:4*N,N/2].reshape( (N*4*N,) )
  iz_slice=iz[0:N,0:4*N,N/2].reshape( (N*4*N,) )
  
  print tau
  
  pyplot.figure(figsize=(12,12))
  
  instance.evolve(0.25*tau) 
  print instance.get_time()
  
  rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(ix_slice,iy_slice,iz_slice)
  pyplot.subplot(4,1,1)
  pyplot.imshow( rho.reshape((N,4*N))**0.5 )
  
  instance.evolve(1.*tau) 
  print instance.get_time()

  rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(ix_slice,iy_slice,iz_slice)
  pyplot.subplot(4,1,2)
  pyplot.imshow( rho.reshape((N,4*N))**0.5 )

  instance.evolve(1.75*tau)
  print instance.get_time()

  rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(ix_slice,iy_slice,iz_slice)
  pyplot.subplot(4,1,3)
  pyplot.imshow( rho.reshape((N,4*N))**0.5 )

  instance.evolve(2.5*tau)
  print instance.get_time()
  
  rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(ix_slice,iy_slice,iz_slice)
  pyplot.subplot(4,1,4)
  pyplot.imshow( rho.reshape((N,4*N))**0.5 )
  pyplot.savefig('test.png')
