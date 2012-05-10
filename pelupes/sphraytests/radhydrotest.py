import numpy
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.units import constants
from amuse.community.simplex.interface import SimpleX
from amuse.community.fi.interface import Fi
from amuse.community.sphray.interface import SPHRay
from amuse.support.data import core

from amuse.ext.evrard_test import uniform_random_unit_cube,uniform_unit_sphere,glass_unit_cube

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

from amuse.support.io import read_set_from_file
from amuse.support.io import write_set_to_file 

import cPickle

write_snapshots=True

gamma=5./3.
mu=1.| units.amu
muion =0.5 | units.amu
xtrans=0.06
mutrans=mu/(1+xtrans)
Tinit=100. | units.K
Ttrans=13500. | units.K 
Tion=13500. | units.K
rhoinit=0.001 | (units.amu / units.cm**3)
uinit=1/(gamma-1)*constants.kB * Tinit/mu
utrans=1/(gamma-1)*constants.kB * Ttrans/mutrans
uion=1/(gamma-1)*constants.kB * Tion/muion
mp=None

hydro_code=Fi
rad_code=SPHRay

def fake_internal_energy_from_xion(xion):
  u=uinit+(utrans-uinit)*xion/xtrans
  a=numpy.where(xion > xtrans )[0]
  u[a]=utrans+(uion-utrans)*(xion[a]-xtrans)/(1.-xtrans)
  return u

def iliev_test_5_ic( N=10000,
                     Ns=10,
                     L=15. | units.kpc ): 
  print "Initializing iliev_test_5"
  
  mp=rhoinit*(2*L)**3/N

  try:
    f=open("glass%9.9i.pkl"%N,"rb")
    x,y,z=cPickle.load(f)
    f.close()
  except:  
    x,y,z=glass_unit_cube(N,target_rms=0.05).make_xyz()
    f=open("glass%9.9i.pkl"%N,"wb")
    cPickle.dump((x,y,z),f)
    f.close()
    
  p=core.Particles(N)

# set particles homogeneously in space
  p.x=L*x
  p.y=L*y
  p.z=L*z

# set other properties
  p.h_smooth=0. | units.parsec
  p.vx= 0. | (units.km/units.s)
  p.vy= 0. | (units.km/units.s)
  p.vz= 0. | (units.km/units.s)
  p.u= uinit
  p.rho = rhoinit
  p.mass=mp
  p.flux=0. | (units.s**-1)
  p.xion=0. | units.none

  sources=core.Particles(Ns)
  x,y,z=uniform_unit_sphere(Ns).make_xyz()
  if Ns == 1:
    x,y,z=0.,0.,0.

  sources.x=L*x*(1./N)**(1./3)/10
  sources.y=L*y*(1./N)**(1./3)/10
  sources.z=L*z*(1./N)**(1./3)/10
  sources.luminosity=(5.e48/Ns) | (units.s**-1)
  sources.SpcType=1.

  return p,sources

def iliev_test_5( N=10000,
                  Ns=10,
                  L=15. | units.kpc,
                  dt=1|units.Myr, rad_parameters=dict()):
  
  rad_parameters["box_size"]=2*L
  
  gas,sources=iliev_test_5_ic(N,Ns,L)                

  conv=nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
     
  sph=Fi(conv,mode='periodic')   

  sph.parameters.use_hydro_flag=True
  sph.parameters.radiation_flag=False
  sph.parameters.self_gravity_flag=False
  sph.parameters.gamma=1
  sph.parameters.isothermal_flag=True
  sph.parameters.integrate_entropy_flag=False
  sph.parameters.timestep=dt/2  
  sph.parameters.verbosity=0
  sph.parameters.periodic_box_size=2*L
  sph.gas_particles.add_particles(gas)
  print sph.parameters.timestep.in_(units.Myr)

  rad=SPHRay()#redirection='none')

  for x in rad_parameters:
    print x, rad_parameters[x]
    rad.parameters.__setattr__(x,rad_parameters[x])

#  gas.u=100*gas.u
  rad.gas_particles.add_particles(gas)
  rad.src_particles.add_particles(sources)

                                    
  return sph,rad

def update_pos_rho(sys,pa):
  p=pa.copy()
  p.rho=3./4/numpy.pi*8*p.mass/p.h_smooth**3
  channel = p.new_channel_to(sys.gas_particles)
  channel.copy_attributes(["x","y","z","rho","h_smooth"])

def update_u(sys,pa):
  sys.particles.u=pa.u

def update_u_dudt(sys,pa):
  sys.particles.u=pa.u
#  sys.particles.du_dt=pa.du_dt
  
def update_u_fake(sys,pa):
  p=pa.copy()
  p.u=fake_internal_energy_from_xion(p.xion)
  channel = p.new_channel_to(sys.particles)
  channel.copy_attributes(["u"])  

def rad_only_evolve(sph,rad,tend,dt):
  global label

  t=rad.model_time
  while t<tend-dt/2:
    print t.in_(units.Myr)        
    rad.evolve_model(t+dt)
    
    t+=dt
    i=t.value_in(units.Myr)

    if i%10==0: 
      p=sph.particles.copy()
      channel = rad.gas_particles.new_channel_to(p)
      channel.copy_attributes(["xion","rho","u"])  
      print (2./3.*mu/constants.kB*p.u.amin()).in_(units.K), (2./3.*mu/constants.kB*p.u.amax()).in_(units.K)
      if write_snapshots:
        write_set_to_file(p,label+"-%6.6i"%int(i),"amuse",
                          append_to_file=False)    
  return t

def radhydro_evolve_fake(sph,rad,tend,dt):
  global label

  t=rad.model_time  
  while t<tend-dt/2:
    print t.in_(units.Myr)        
    print "sph1"
    sph.evolve_model(t+dt/2)    
    print "rad"
    update_pos_rho(rad,sph.gas_particles)
    rad.evolve_model(t+dt)
    update_u_fake(sph,rad.gas_particles)
    print (2./3.*mu/constants.kB*sph.gas_particles.u.amin()).in_(units.K), (2./3.*mu/constants.kB*sph.gas_particles.u.amax()).in_(units.K)
    print "sph2"
    sph.evolve_model(t+dt)    
    
    t+=dt
    i=t/dt
    if int(i)%10==0:
      if write_snapshots:
        p=rad.gas_particles.copy()
        channel = sph.gas_particles.new_channel_to(p)
        channel.copy_attributes(["x","y","z","vx","vy","vz","rho","u"])  	  
        write_set_to_file(p,label+"-%6.6i"%int(i),"amuse",
                        append_to_file=False)    
  return t

def radhydro_evolve(sph,rad,tend,dt):
  global label

  t=rad.model_time  
  while t<tend-dt/2:
    print t.in_(units.Myr)        
    print "sph1"
    sph.evolve_model(t+dt/2)    
    print "rad"
    update_pos_rho(rad,sph.gas_particles)
    rad.evolve_model(t+dt)
    update_u(sph,rad.gas_particles)
    print (2./3.*mu/constants.kB*sph.gas_particles.u.amin()).in_(units.K), (2./3.*mu/constants.kB*sph.gas_particles.u.amax()).in_(units.K)
    print "sph2"
    sph.evolve_model(t+dt)    
    
    t+=dt
    i=t/dt
    if int(i)%10==0:
      if write_snapshots:
        p=rad.gas_particles.copy()
        channel = sph.gas_particles.new_channel_to(p)
        channel.copy_attributes(["x","y","z","vx","vy","vz","rho","u"])  	  
        write_set_to_file(p,label+"-%6.6i"%int(i),"amuse",
                        append_to_file=False)    
  return t

def average(N,r,dat):
  n=len(r)
  a=numpy.argsort(r)
  i=0
  r_a=[]
  dat_a=[]
  while i < n:
    ra=r[a[i:i+N]].sum()/min(n-i,N)
    da=dat[a[i:i+N]].sum()/min(n-i,N)
    r_a.append(ra)
    dat_a.append(da)
    i=i+N
  return numpy.array(r_a),numpy.array(dat_a)

def aplot(i, tag, xyc,xlim=None,ylim=None,ylabel=""):
  pyplot.figure(figsize=(8,8))
  for x,y,c in xyc:
    xa,ya=average(100,x,y)
    pyplot.semilogy(xa,ya,c)  
  pyplot.xlabel('L/Lbox')  
  pyplot.ylabel(ylabel)  
  if xlim is not None:
    pyplot.xlim(xlim)
  if ylim is not None:
    pyplot.ylim(ylim)
  pyplot.savefig(tag+'-%6.6i.png'%i)


def test(evolver, rad_parameters=dict(),dt=1. | units.Myr ):    
  N=100000
  Ns=1
  L=15. | units.kpc
  tplot=[10.| units.Myr,200.| units.Myr,500.| units.Myr]

  rad_parameters["box_size"]=2*L

  sph,rad=iliev_test_5(N,Ns,L,dt,rad_parameters=rad_parameters)

  p=sph.particles.copy()
  channel = rad.gas_particles.new_channel_to(p)
  channel.copy_attributes(["xion","rho","u"])  
  print (2./3.*mu/constants.kB*p.u.amin()).in_(units.K), (2./3.*mu/constants.kB*p.u.amax()).in_(units.K)
  if write_snapshots:
    write_set_to_file(p,label+"-%6.6i"%0,"amuse",
                      append_to_file=False)    

  for tend in tplot:
    evolver(sph,rad,tend,dt)

    p=rad.gas_particles.copy()
    channel = sph.gas_particles.new_channel_to(p)
    if evolver==rad_only_evolve:
      channel.copy_attributes(["vx","vy","vz"])  
    else:
      channel.copy_attributes(["vx","vy","vz","u"])  

    r=(((p.x)**2+(p.y)**2+(p.z)**2) ** 0.5).value_in(units.kpc)
    v=((p.vx**2+p.vy**2+p.vz**2)**0.5).value_in(units.kms)
    cs=(gamma*(gamma-1)*p.u**0.5).value_in(units.kms)
    xion=p.xion
    t=(gamma-1)*mu/(1.+xion)/constants.kB*p.u
    t=t.value_in(units.K)
    dens=(p.rho).value_in(units.amu/units.cm**3)
    pres=(p.u*p.rho).value_in( units.g/units.cm/units.s**2)
    mach=v/cs+1.e-10

    i=int(tend/dt)

    pyplot.figure( figsize=(8,6) ) 
    pyplot.plot( r,xion,'r.')
    pyplot.savefig('rx-%6.6i.png'%i)

    pyplot.figure( figsize=(8,6) ) 
    pyplot.plot( r,t,'r.')
    pyplot.savefig('rt-%6.6i.png'%i)

    aplot(i,label+'-xion',((r/15,xion,'r'),(r/15,1-xion,'g')),
            xlim=(0.,1.),ylim=(1.e-6,1.),ylabel='x,(1-x)')
    aplot(i,label+'-temp',((r/15,t,'r'),),
            xlim=(0.,1.),ylim=(100,1.e5),ylabel="T (K)")
    aplot(i,label+'-pres',((r/15,pres,'r'),),
            xlim=(0.,1.),ylim=(1.e-17,1.e-14),ylabel="pres (g/cm/s**2)")
    aplot(i,label+'-rho',((r/15,dens,'r'),),
            xlim=(0.,1.),ylim=(0.0001,0.01),ylabel='density (amu/cm**3)')
    aplot(i,label+'-mach',((r/15,mach,'r'),),
            xlim=(0.,1.),ylim=(1.e-5,10.),ylabel='Mach')
    aplot(i,label+'-vel',((r/15,v,'r'),),
            xlim=(0.,1.),ylim=(0.1,15.),ylabel='v (km/s)')


if __name__=="__main__":
  global label
  import time
  
  """
  t1=time.time()
  label="test1"
  test(rad_only_evolve)
  t2=time.time()
  print "test1: monochromatic, non dynamic, no thermal solver",
  print t2-t1,"sec"

  t1=time.time()
  label="test2"
  test(rad_only_evolve, dict(number_of_freq_bins=5, blackbody_spectrum_flag=1) )
  t2=time.time()
  print "test2: 5 freq, non dynamic, no thermal solver",
  print t2-t1,"sec"
  
  t1=time.time()
  label="test3"
  test(rad_only_evolve, dict(number_of_freq_bins=5, blackbody_spectrum_flag=1,thermal_evolution_flag=1) )
  t2=time.time()
  print "test3: 5 freq, non dynamic, with thermal solver",
  print t2-t1,"sec"

  
  t1=time.time()
  label="test4"
  test(radhydro_evolve_fake,dict(isothermal_flag=True, number_of_rays=1000,hydrogen_case_A_flag=False))
  t2=time.time()
  print "test4: monochromatic, dynamic, no thermal solver",
  print t2-t1,"sec"
  
  t1=time.time()
  label="test5"
  test(radhydro_evolve_fake, dict(number_of_freq_bins=5, blackbody_spectrum_flag=1) )
  t2=time.time()
  print "test5: 5 freq, dynamic, no thermal solver",
  print t2-t1,"sec"
  
  """
  t1=time.time()
  label="test6"
  test(radhydro_evolve, dict(isothermal_flag=False, number_of_rays=10000,hydrogen_case_A_flag=False) )
  t2=time.time()
  print "test6: 5 freq, dynamic, thermal solver",
  print t2-t1,"sec"
  
