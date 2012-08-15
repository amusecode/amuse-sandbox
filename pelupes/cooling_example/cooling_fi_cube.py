import numpy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

from amuse.units import generic_unit_converter
from amuse.units import units, constants

from amuse.community.fi.interface import Fi

from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.molecular_cloud import ism_cube
from cooling import evolve_internal_energy, global_mu

import cPickle

#import logging
#logging.basicConfig(level=logging.DEBUG)

def make_map(sph,L,N=100):

    x,y=numpy.indices( ( N+1,N+1 ))

    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=x*0.
        
    vx=units.kms(numpy.zeros_like(x.number))
    vy=units.kms(numpy.zeros_like(x.number))
    vz=units.kms(numpy.zeros_like(x.number))

    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    kin_energy=0.5*(rhovx**2+rhovy**2+rhovz**2)/rho
    internal_energy=abs(rhoe-kin_energy)/rho
    temp=global_mu / constants.kB * internal_energy

    rho=rho.reshape((N+1,N+1))
    temp=temp.reshape((N+1,N+1))

    return numpy.transpose(rho),numpy.transpose(temp)

if __name__=="__main__":
  import time

  mapstep=1
  N=10000
  L= 10| units.parsec
  rho= 1.14*1. | units.amu/units.cm**3
  u=(5.e11 | units.erg/units.g).in_(units.cm**2/units.s**2) # =5000K
  
  tend=0.6 | units.Myr
  dt=2500/mapstep | units.yr
  
  print (u**0.5).in_(units.kms)
  print ((L/u**0.5)/dt)
    
  particles=ism_cube(N, L,rho,u,seed=123456,base_grid=regular_grid_unit_cube,nf=64).result

  print len(particles)

  UnitLength=L
  UnitMass=particles.mass.sum()
  UnitVelocity=units.kms
  
  convert=generic_unit_converter.ConvertBetweenGenericAndSiUnits(UnitLength, UnitMass, UnitVelocity)
  
  sph=Fi(convert,mode='periodic',redirection="none")#number_of_workers=1,use_gl=False,debugger='gdb')
  
  sph.parameters.periodic_box_size=2*L
  sph.parameters.use_hydro_flag=True
  sph.parameters.radiation_flag=False
  sph.parameters.self_gravity_flag=False
  sph.parameters.gamma=1
  sph.parameters.isothermal_flag=True
  sph.parameters.integrate_entropy_flag=False
  sph.parameters.timestep=dt/2  
  sph.parameters.verbosity=0
  
  sph.gas_particles.add_particles(particles)
  
  label="test-iso"
  
#  sph.start_viewer()
  
  i=0
  step=0
  t=0. | units.Myr
  t1=time.time()
  
  times=[]
  densfrac=[]
  coolfrac=[]
  meant=[]
  
  while t<(tend-dt/2):
    step=step+1
    sph.evolve_model(t + dt/2.0)
    print (global_mu / constants.kB * sph.gas_particles.u.amin()).in_(units.K), 
    print (global_mu / constants.kB * sph.gas_particles.u.mean()).in_(units.K), 
    print (global_mu / constants.kB * sph.gas_particles.u.amax()).in_(units.K)
    sph.gas_particles.u = evolve_internal_energy(sph.gas_particles.u, dt, 
        sph.gas_particles.rho/global_mu, sph.gas_particles.du_dt)
    print sph.gas_particles.rho.amin().in_(units.amu/units.cm**3),
    print sph.gas_particles.rho.amax().in_(units.amu/units.cm**3)
    print (global_mu / constants.kB * sph.gas_particles.u.amin()).in_(units.K), 
    print (global_mu / constants.kB * sph.gas_particles.u.mean()).in_(units.K), 
    print (global_mu / constants.kB * sph.gas_particles.u.amax()).in_(units.K)
    sph.evolve_model(t + dt)
    t=t+dt
    print t.in_(units.Myr),sph.model_time.in_(units.Myr)

    if step%mapstep == 0:
      i=i+1
      rho,temp=make_map(sph,2*L,N=100)
      
      f=pyplot.figure(figsize=(8,8))
      LL=L.number
      pyplot.imshow(numpy.log10(temp.value_in(units.K)),
          extent=[-LL,LL,-LL,LL],vmin=2,vmax=5,origin='lower')
      pyplot.xlabel("parsec")
      pyplot.savefig('tmap_fi-'+label+'-%6.6i.png'%i)
      f.clear()
      pyplot.close(f)
      print 'tmap-%6.6i.png'%i
      
      f=pyplot.figure(figsize=(8,8))
      LL=L.number
      pyplot.imshow(numpy.log10(rho.value_in(units.amu/units.cm**3)),
          extent=[-LL,LL,-LL,LL],vmin=-1,vmax=2,origin='lower')
      pyplot.xlabel("parsec")
      pyplot.savefig('map_fi-'+label+'-%6.6i.png'%i)
      f.clear()
      pyplot.close(f)
      print 'map-%6.6i.png'%i
      
      times.append(t.value_in(units.Myr))
      rho=sph.particles.rho
      a=numpy.where(rho > (2. | units.amu/units.cm**3))[0]
      if len(a) > 0:
        frac=len(a)/(1.*len(rho))
      else:
        frac=0.
      densfrac.append(frac)
      
      temp=(global_mu / constants.kB * sph.gas_particles.u)
      mt=temp.mean()
      a=numpy.where(temp < (1000. | units.K))[0]
      if len(a) > 0:
        frac=len(a)/(1.*len(temp))
      else:
        frac=0.
      coolfrac.append(frac)
      meant.append(mt.value_in(units.K))
    
  f=pyplot.figure(figsize=(8,6))
  pyplot.plot(times,densfrac,'r')
  pyplot.plot(times,coolfrac,'g')
  pyplot.xlabel("time (Myr)")
  pyplot.ylabel("dens, cool mass fraction")
  pyplot.ylim(0,1)
  pyplot.xlim(0,0.6)
  pyplot.savefig('denscool-'+label+'.eps')
  f.clear()
  pyplot.close(f)

  f=pyplot.figure(figsize=(8,6))
  pyplot.plot(times,meant,'g')
  pyplot.xlabel("time (Myr)")
  pyplot.ylabel("T (K)")
  pyplot.ylim(0,10000.)
  pyplot.xlim(0,0.6)
  pyplot.savefig('meant-'+label+'.eps')
  f.clear()
  pyplot.close(f)
  
  f=open("fractions-"+label+".pkl","wb")
  cPickle.dump([times,densfrac,coolfrac,meant],f)
  f.close()

  t2=time.time()
  print t2-t1
    
