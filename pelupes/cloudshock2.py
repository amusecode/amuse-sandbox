import os
import sys

from amuse.legacy.capreole.interface import Capreole
from amuse.legacy.capreole.interface import GLCapreole

import numpy
from matplotlib import pyplot
import matplotlib.cm as cm

def cloud_shock_problem(N=32,NN=4,M=10.,xi=10.,L=10.,rc=1.,gamma=5./3):

  cs=numpy.sqrt(gamma*(gamma-1))
  cs_out=numpy.sqrt(gamma*(gamma-1)*xi)
  vs=cs_out*M

  dl=L/N
  dr=numpy.sqrt(3.)*dl

  cx=dl*0.4
  cy=dl*0.
  cz=dl*0.6

  ix,iy,iz=numpy.indices( (N,4*N,N) )

  ix=ix.reshape((4*N*N*N,))+1
  iy=iy.reshape((4*N*N*N,))+1
  iz=iz.reshape((4*N*N*N,))+1

  x=L*(ix/(1.*N)-0.5-1/2./N)
  y=L*(iy/(1.*N)-0.5-1/2./N)
  z=L*(iz/(1.*N)-0.5-1/2./N)

  r=((x-cx)**2+(y-cy)**2+(z-cz)**2)**0.5

  rho=numpy.zeros_like(r)
  rhvx=numpy.zeros_like(r)
  rhvy=numpy.zeros_like(r)
  rhvz=numpy.zeros_like(r)
  en=numpy.zeros_like(r)

  selection= r > (rc+dr)
  selection=numpy.compress(selection,numpy.arange(len(r))) 
  rho[selection]=1./xi
  rhvx[selection]=0.
  rhvy[selection]=1./xi*vs
  rhvz[selection]=0.
  en[selection]=xi+1./xi*vs**2
  
  selection= r < (rc-dr)
  selection=numpy.compress(selection,numpy.arange(len(r))) 
  rho[selection]=1.
  rhvx[selection]=0.
  rhvy[selection]=0.
  rhvz[selection]=0.
  en[selection]=1.
    
  ix,iy,iz=numpy.indices( (NN,NN,NN) )
  ix=ix.reshape((NN*NN*NN,))+1
  iy=iy.reshape((NN*NN*NN,))+1
  iz=iz.reshape((NN*NN*NN,))+1

  dx=dl*(ix/(1.*NN)-0.5-1/2./NN)
  dy=dl*(iy/(1.*NN)-0.5-1/2./NN)
  dz=dl*(iz/(1.*NN)-0.5-1/2./NN)

  drho=numpy.zeros_like(dx)
  drhvx=numpy.zeros_like(dx)
  drhvy=numpy.zeros_like(dx)
  drhvz=numpy.zeros_like(dx)
  den=numpy.zeros_like(dx)

  selection= numpy.logical_and( r>=(rc-dr), r <= (rc+dr))
  selection=numpy.compress(selection,numpy.arange(len(r))) 

  for i in selection:
    xx=x[i]  
    yy=y[i]
    zz=z[i]
    r=((xx+dx-cx)**2+(yy+dy-cy)**2+(zz+dz-cz)**2)**0.5
    dselection=r>rc
    dselection=numpy.compress(dselection,numpy.arange(len(r))) 
    drho[dselection]=1./xi
    drhvx[dselection]=0.
    drhvy[dselection]=1./xi*vs
    drhvz[dselection]=0.
    den[dselection]=xi+1./xi*vs**2
    
    dselection=r<=rc
    dselection=numpy.compress(dselection,numpy.arange(len(r))) 
    drho[dselection]=1.
    drhvx[dselection]=0.
    drhvy[dselection]=0.
    drhvz[dselection]=0.
    den[dselection]=1.

    rho[i]=drho.mean()
    rhvx[i]=drhvx.mean()
    rhvy[i]=drhvy.mean()
    rhvz[i]=drhvz.mean()
    en[i]=den.mean()
    
  rho=rho.reshape((N,4*N,N))
  rhvx=rhvx.reshape((N,4*N,N))
  rhvy=rhvy.reshape((N,4*N,N))
  rhvz=rhvz.reshape((N,4*N,N))
  en=en.reshape((N,4*N,N))

  ix,iy,iz=numpy.indices( (N,4*N,N) )

  ix=ix+1
  iy=iy+1
  iz=iz+1

  return ix,iy,iz,rho,rhvx,rhvy,rhvz,en

def cloudshock(N=32,NN=4,xi=10.,M=2.7,L=10.,rc=1.,gamma=5./3):
  ix,iy,iz,rho,rhvx,rhvy,rhvz,en=cloud_shock_problem(N=N,NN=NN,M=M,xi=xi,L=L,rc=rc,gamma=gamma)
  ix=ix.reshape( (4*N*N*N,) )
  iy=iy.reshape( (4*N*N*N,) )
  iz=iz.reshape( (4*N*N*N,) )
  rho=rho.reshape( (4*N*N*N,) )
  rhvx=rhvx.reshape( (4*N*N*N,) )
  rhvy=rhvy.reshape( (4*N*N*N,) )
  rhvz=rhvz.reshape( (4*N*N*N,) )
  en=en.reshape( (4*N*N*N,) )
  instance=Capreole(number_of_workers=3)
#  instance=Capreole(name_of_the_worker="worker")
  instance.initialize_code()
  instance.set_boundary("periodic","periodic","interface","outflow",
                          "periodic","periodic")
  instance.set_boundary_innerystate(rho[0],rhvx[0],rhvy[0],rhvz[0],en[0])                        
  instance.setup_mesh(N,4*N,N,L,4*L,L)
  instance.commit_parameters()
  instance.set_grid_state(ix,iy,iz,rho,rhvx,rhvy,rhvz,en)
  instance.initialize_grid(0.0)
  if hasattr(instance,'viewer'):
    instance.viewer()
  return ix,iy,iz,instance

if __name__=="__main__":
  N=160
  
  rc=1.
  M=2.7
  xi=10.
  gamma=5./3

  ix,iy,iz,instance=cloudshock(N=N,rc=rc,M=M,xi=xi)
  
  cs=numpy.sqrt(gamma*(gamma-1))
  cs_out=numpy.sqrt(gamma*(gamma-1)*xi)
  vs=cs_out*M
  tau=1.6*2*rc*xi**0.5/vs
  
  ix=ix.reshape( (N,4*N,N) )
  iy=iy.reshape( (N,4*N,N) )
  iz=iz.reshape( (N,4*N,N) )
  
  ix_slice=ix[0:N,0:4*N,N/2].reshape( (N*4*N,) )
  iy_slice=iy[0:N,0:4*N,N/2].reshape( (N*4*N,) )
  iz_slice=iz[0:N,0:4*N,N/2].reshape( (N*4*N,) )
  
  print tau
  

  scl=numpy.log
  pyplot.figure(figsize=(12,12))
  
  instance.evolve(0.25*tau) 
  
  
  rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(ix_slice,iy_slice,iz_slice)
  pyplot.subplot(4,1,1)
  pyplot.imshow( scl(rho.reshape((N,4*N))),vmin=scl(1./xi))
  
  instance.evolve(0.75*tau) 

  rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(ix_slice,iy_slice,iz_slice)
  pyplot.subplot(4,1,2)
  pyplot.imshow( scl(rho.reshape((N,4*N))),vmin=scl(1./xi))

  instance.evolve(1.5*tau)

  rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(ix_slice,iy_slice,iz_slice)
  pyplot.subplot(4,1,3)
  pyplot.imshow( scl(rho.reshape((N,4*N))),vmin=scl(1./xi))

  instance.evolve(2.25*tau)
  
  rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(ix_slice,iy_slice,iz_slice)
  pyplot.subplot(4,1,4)
  pyplot.imshow( scl(rho.reshape((N,4*N))),vmin=scl(1./xi))
  pyplot.savefig('test.png')

