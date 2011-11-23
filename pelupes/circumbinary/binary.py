import os
import sys
import numpy
import cPickle
import time

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

from amuse.community.twobody.twobody import TwoBody

from amuse.units import units
from amuse.units import nbody_system
from amuse.units import constants

from amuse.datamodel import Particles

#from amuse.community.fi.interface import Fi
from boxedfi import BoxedFi as Fi

from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.ext.evrard_test import regular_grid_unit_cube

from fast import FAST

from amuse.io import write_set_to_file

from directsum import directsum

def binary(interface,m1=1.|units.MSun,m2=.0001| units.MSun,r1=None,r2=None,ecc=0,P=1| units.yr):

  mu=constants.G*(m1+m2)
  a=(P/(2*numpy.pi)*mu**0.5)**(2./3.)

  f1=m2/(m1+m2)
  f2=m1/(m1+m2)

  rmax=a*(1+ecc)

  r0=rmax
  print a.in_(units.AU),r0.in_(units.AU)

  
  h=(a*mu*(1-ecc**2))**0.5
  v0=h/r0

  bin=Particles(2)

  bin[0].mass=m1
  bin[0].x=r0*f1
  bin[0].vy=v0*f1
  bin[1].mass=m2
  bin[1].x=-r0*f2
  bin[1].vy=-v0*f2

  bin.y=0*r0
  bin.z=0.*r0
  bin.vx=0*v0
  bin.vz=0.*v0
  if r1 is None:
    bin[0].radius=(1.|units.RSun)*(m1/(1.|units.MSun))**(1./3.)
  else:
    bin[0].radius=r1
  if r2 is None:
    bin[1].radius=(1.|units.RSun)*(m2/(1.|units.MSun))**(1./3.)
  else:
    bin[1].radius=r2
  radius=[0.,0.]

  convert=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)

  nb = interface(convert,redirection="none")
  
  nb.particles.add_particles(bin)

  return nb

def sphdisc(interface,N=10000,Mstar=1| units.MSun, Rmin=1|units.AU,Rmax=10|units.AU):

    convert=nbody_system.nbody_to_si(Mstar, Rmin.unit)
    proto=ProtoPlanetaryDisk(N,convert_nbody=convert,densitypower=1.5,
                               Rmin=Rmin.number,Rmax=Rmax.number,q_out=1.5)
    gas=proto.result
    gas.h_smooth=0.06 | units.AU
    gas.u0=gas.u.copy()

    sph=interface(convert,redirection='none')

    sph.parameters.use_hydro_flag=True
    sph.parameters.radiation_flag=False
    sph.parameters.self_gravity_flag=True
    sph.parameters.gamma=1.
    sph.parameters.isothermal_flag=True
    sph.parameters.integrate_entropy_flag=False
    sph.parameters.timestep=1. | units.day  
    sph.parameters.verbosity=0
    sph.parameters.courant=0.2    

#    print sph.parameters
#    print sph.parameters.periodic_box_size.in_(units.AU)
    print gas.mass.sum().in_(units.MSun)

    sph.gas_particles.add_particles(gas)

    return sph,gas

def make_map(sph,N=100,L=1):

    x,y=numpy.indices( ( N+1,N+1 ))

    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=x*0.
    vx=0.*x
    vy=0.*x
    vz=0.*x

    x=units.AU(x)
    y=units.AU(y)
    z=units.AU(z)
    vx=units.kms(vx)
    vy=units.kms(vy)
    vz=units.kms(vz)

    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    rho=rho.reshape((N+1,N+1))

    return numpy.transpose(rho)

def make_phi_map(sph,N=100,Rrange=(0.3,2),phioffset=0.):


    phi,r=numpy.mgrid[0:2*numpy.pi:N*1j,Rrange[0]:Rrange[1]:N*1j] 
    phi=phi.flatten()
    r=r.flatten()

    x,y=r*numpy.cos(phi+phioffset),r*numpy.sin(phi+phioffset)
    z=x*0.

    vx=0.*x
    vy=0.*x
    vz=0.*x

    x=units.AU(x)
    y=units.AU(y)
    z=units.AU(z)
    vx=units.kms(vx)
    vy=units.kms(vy)
    vz=units.kms(vz)

    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    rho=rho.reshape((N,N))

    return numpy.transpose(rho)

def handle_eos(gas,memgas,rhotrans=(1.e-5 | units.g/units.cm**3), gamma=1.4):
  channel=gas.new_channel_to(memgas)
  channel.copy_attribute("rho")
  a=memgas.select_array(lambda rho: rho> rhotrans,["rho"])
  a.u=a.u0*(a.rho/rhotrans)**gamma
  a=memgas.select_array(lambda rho: rho<= rhotrans,["rho"])
  a.u=a.u0
  channel=memgas.new_channel_to(gas)
  channel.copy_attribute("u")

def sink_particles(sinks,sources,Raccretion=0.1 | units.AU):

  R2=Raccretion**2
  for s in sinks:
     ms=s.mass
     xs,ys,zs=s.x,s.y,s.z
     insink=sources.select_array(lambda x,y,z: (x-xs)**2+(y-ys)**2+(z-zs)**2 < R2,['x','y','z'])  
     if len(insink) > 0:
       cm=s.position*ms
       p=s.velocity*ms

       s.mass+=insink.total_mass()
       s.position=(cm+insink.center_of_mass()*insink.total_mass())/s.mass
       s.velocity=(p+insink.total_momentum())/s.mass
# we lose angular momentum !    
       sources.remove_particles(insink)   



if __name__=="__main__":


  Raccretion=0.1 | units.AU
  pplanet=228.776 #  | units.day
  bin=binary(TwoBody,m1=0.6897 |units.MSun,m2=.20255| units.MSun,
    r1=0.6489 | units.RSun,r2=0.22623 | units.RSun,
    ecc=0.15944,P=41.08| units.day)
    
  disc,gas=sphdisc(Fi,100000,bin.particles.mass.sum(),Rmin=0.5|units.AU,Rmax=8|units.AU)
  
  dt_int= 2. | units.day
  
  directsum_disc=directsum( (disc,) )
    
  bridge=FAST(verbose=False)
  bridge.set_timestep(dt_int)
  bridge.add_system(bin, (directsum_disc,), False)
  bridge.add_system(disc, (bin,), False)

  tnow=0. |  units.day  
  dt=4. | units.day
  tend=25. | units.yr

  L=6.

  i=0
  while tnow < tend-dt/2:

    if i%20==0:
      write_set_to_file(bin.particles,'snap/bin-%6.6i'%(i),'amuse')
      write_set_to_file(disc.gas_particles,'snap/disc-%6.6i'%(i),'amuse')

    tnow+=dt
    i+=1
    handle_eos(disc.particles,gas,rhotrans=(1.e-5 | units.g/units.cm**3))
    sink_particles(bin.particles,disc.particles,Raccretion=Raccretion)
    bridge.evolve_model(tnow)
    print tnow

    if i%1==0:
      rho=make_map(disc,N=200,L=L)
      f=pyplot.figure(figsize=(8,8))
      pyplot.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)),
          extent=[-L/2,L/2,-L/2,L/2],vmin=12,vmax=17.5,origin='lower')    
      pyplot.plot(bin.particles.x.value_in(units.AU),
                  bin.particles.y.value_in(units.AU),'r+')
      pyplot.xlim(-L/2,L/2)
      pyplot.ylim(-L/2,L/2)
      pyplot.title(tnow)
      pyplot.xlabel('AU')
      pyplot.savefig('map/map%6.6i.png'%i)
      f.clear()
      pyplot.close(f)

      offset=numpy.mod(2*numpy.pi*tnow.value_in(units.day)/pplanet,2*numpy.pi)

      rho=make_phi_map(disc,N=200,phioffset=offset)
      f=pyplot.figure(figsize=(12,4))
      pyplot.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)),
          extent=[0,2*numpy.pi,0.3,2],vmin=12,vmax=17.5,origin='lower')    
      x=bin.particles.x.value_in(units.AU)
      y=bin.particles.y.value_in(units.AU)
      r=(x**2+y**2)**0.5
      phi=numpy.arctan2(y,x)-offset
      phi=numpy.mod(phi,2*numpy.pi)
      pyplot.plot(phi,r,'r+')
      pyplot.xlabel('phi')
      pyplot.ylabel('R')
      pyplot.xlim(0,2.*numpy.pi)
      pyplot.ylim(0.,2)
      pyplot.savefig('map/phi-%6.6i.png'%i)
      f.clear()
      pyplot.close(f)
           
