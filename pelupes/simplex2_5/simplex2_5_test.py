import numpy
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.units import constants
from amuse.community.simplex2_5.interface import SimpleX
from amuse.community.fi.interface import Fi
from amuse.support.data import core

from amuse.ext.evrard_test import uniform_random_unit_cube,uniform_unit_sphere

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

from amuse.support.io import read_set_from_file
from amuse.support.io import write_set_to_file 

write_snapshots=False


mu=1.| units.amu
muion =0.5 | units.amu
xtrans=0.06
mutrans=mu/(1+xtrans)
Tinit=100. | units.K
Ttrans=13500. | units.K 
Tion=13500. | units.K
rhoinit=0.001 | (units.amu / units.cm**3)
uinit=constants.kB * Tinit/mu
utrans=constants.kB * Ttrans/mutrans
uion=constants.kB * Tion/muion
mp=None

def fake_internal_energy_from_xion(xion):
#  return uinit+(uion-uinit)*xion 
  u=uinit+(utrans-uinit)*xion/xtrans
  a=numpy.where(xion > xtrans )[0]
  u[a]=utrans+(uion-utrans)*(xion[a]-xtrans)/(1.-xtrans)
  return u
  
def glass(N, target_rms=0.01):
  
  try:
    particles=read_set_from_file("glass%9.9i.dump"%N,"amuse")
  except:
    if target_rms < 0.001:
      print "warning: target_rms highly unlikely to succeed"

    L=1| nbody_system.length
    dt=0.01 | nbody_system.time
    x,y,z=uniform_random_unit_cube(N).make_xyz()
    vx,vy,vz=uniform_unit_sphere(N).make_xyz()


    p=core.Particles(N)
    p.x=L*x
    p.y=L*y
    p.z=L*z
    p.h_smooth=0. * L
    p.vx= 0.1*vx | (nbody_system.speed)
    p.vy= 0.1*vy | (nbody_system.speed)
    p.vz= 0.1*vz | (nbody_system.speed)
    p.u= (0.1*0.1) | nbody_system.speed**2 
    p.mass=(8./N) | nbody_system.mass

    sph=Fi(use_gl=False,mode='periodic',redirection='none')   
    sph.initialize_code()

    sph.parameters.use_hydro_flag=True
    sph.parameters.radiation_flag=False
    sph.parameters.self_gravity_flag=False
    sph.parameters.gamma=1
    sph.parameters.isothermal_flag=True
    sph.parameters.integrate_entropy_flag=False
    sph.parameters.timestep=dt  
    sph.parameters.verbosity=0
    sph.parameters.periodic_box_size=2*L
    sph.parameters.artificial_viscosity_alpha=1.| units.none
    sph.parameters.beta=2. | units.none
    sph.commit_parameters()
    sph.gas_particles.add_particles(p)
    sph.commit_particles()

  #  sph.start_viewer()

    t=0. | nbody_system.time
    rms=1.
    i=0
    while rms > target_rms:
      t=t+(0.25 | nbody_system.time)
      sph.evolve_model(t)
      h=sph.particles.h_smooth.value_in(nbody_system.length)
      rho=h**(-3.)
      rms=rho.std()/rho.mean()
      print rms



    write_set_to_file(sph.particles,"glass%9.9i.dump"%N,"amuse",append_to_file=False)    

    particles=sph.particles.copy()
    del sph
  
  x=particles.x.value_in(nbody_system.length)
  y=particles.y.value_in(nbody_system.length)
  z=particles.z.value_in(nbody_system.length)

  return x,y,z


def iliev_test_5_ic( N=10000,
                  Ns=10,
                  L=15. | units.kpc ):

  mp=rhoinit*(2*L)**3/N
 
  print "Initializing iliev_test_5"

#  x,y,z=uniform_random_unit_cube(N).make_xyz()
  x,y,z=glass(N,target_rms=0.05)
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

  sources.x=L*x*(1./N)**(1./3)/10
  sources.y=L*y*(1./N)**(1./3)/10
  sources.z=L*z*(1./N)**(1./3)/10
  sources.rho=rhoinit/100.
  sources.flux=(5.e48/Ns) | (units.s**-1)
  sources.xion=1. | units.none
  
  return p,sources

def iliev_test_5( N=10000,
                  Ns=10,
                  L=15. | units.kpc,
                  dt=None, rad_parameters=None):
  
  if rad_parameters is None:
    rad_parameters=dict(box_size=2*L)
  
  gas,sources=iliev_test_5_ic(N,Ns,L)                

  conv=nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
     
  sph=Fi(conv,use_gl=False,mode='periodic')   
  sph.initialize_code()

  sph.parameters.use_hydro_flag=True
  sph.parameters.radiation_flag=False
  sph.parameters.self_gravity_flag=False
  sph.parameters.gamma=1
  sph.parameters.isothermal_flag=True
  sph.parameters.integrate_entropy_flag=False
  sph.parameters.timestep=dt  
  sph.parameters.verbosity=0
  sph.parameters.periodic_box_size=2*L
  sph.commit_parameters()
  sph.gas_particles.add_particles(gas)
  sph.commit_particles()

  print sph.parameters.timestep.in_(units.Myr)

#  sph.start_viewer()
         
  rad=SimpleX(number_of_workers = 2,redirection='none')

  for x in rad_parameters:
    print x, rad_parameters[x]
    rad.parameters.__setattr__(x,rad_parameters[x])

  rad.commit_parameters()

  gas.add_particles(sources)
  rad.particles.add_particles(gas)
  rad.commit_particles()
                                    
  return sph,rad


def update_pos_rho(sys,pa):
  p=pa.copy()
  p.rho=3./4/numpy.pi*8*p.mass/p.radius**3
  channel = p.new_channel_to(sys.particles)
  channel.copy_attributes(["x","y","z","rho"])

def update_u(sys,pa):
  sys.particles.u=pa.u
  
def update_u_fake(sys,pa):
  p=pa.copy()
  p.u=fake_internal_energy_from_xion(p.xion)
  channel = p.new_channel_to(sys.particles)
  channel.copy_attributes(["u"])  
  
def radhydro_evolve(sph,rad,tend,dt):
  global label
  i=0
  if write_snapshots:
    write_set_to_file(rad.particles,label+"-%6.6i"%i,"amuse",
                        append_to_file=False)    

  print (mu/constants.kB*sph.gas_particles.u.amin()).in_(units.K), (mu/constants.kB*sph.gas_particles.u.amax()).in_(units.K)
  #raise Exception

  t=0. | units.Myr
  while t<tend-dt/2:
    print t        
    print "sph1"
    sph.evolve_model(t+dt/2)    
    print "rad"
    update_pos_rho(rad,sph.gas_particles)
    rad.evolve_model(t+dt)
    update_u(sph,rad.particles)
    print (mu/constants.kB*sph.gas_particles.u.amin()).in_(units.K), (mu/constants.kB*sph.gas_particles.u.amax()).in_(units.K)
    print "sph2"
    sph.evolve_model(t+dt)    
    
    t+=dt
    i+=1
    
    if write_snapshots:
      write_set_to_file(rad.particles,label+"-%6.6i"%i,"amuse",
                        append_to_file=False)    

def radhydro_evolve_fake(sph,rad,tend,dt):
  global label
  i=0
  if write_snapshots:
    write_set_to_file(rad.particles,label+"-%6.6i"%i,"amuse",
                        append_to_file=False)    

  print (mu/constants.kB*sph.gas_particles.u.amin()).in_(units.K), (mu/constants.kB*sph.gas_particles.u.amax()).in_(units.K)

  t=0. | units.Myr
  while t<tend-dt/2:
    print t        
    print "sph1"
    sph.evolve_model(t+dt/2)    
    print "rad"
    update_pos_rho(rad,sph.gas_particles)
    rad.evolve_model(t+dt)
    update_u_fake(sph,rad.particles)
    print (mu/constants.kB*sph.gas_particles.u.amin()).in_(units.K), (mu/constants.kB*sph.gas_particles.u.amax()).in_(units.K)
    print "sph2"
    sph.evolve_model(t+dt)    
    
    t+=dt
    i+=1
    if write_snapshots:
      write_set_to_file(rad.particles,label+"-%6.6i"%i,"amuse",
                        append_to_file=False)    

 

def rad_only_evolve(sph,rad,tend,dt):
  global label

  i=0
  p=sph.particles.copy()
  channel = rad.particles.new_channel_to(p)
  channel.copy_attributes(["xion","rho","u"])  
  print (mu/constants.kB*p.u.amin()).in_(units.K), (mu/constants.kB*p.u.amax()).in_(units.K)

  if write_snapshots:
    write_set_to_file(p,label+"-%6.6i"%i,"amuse",
                        append_to_file=False)    

  t=0. | units.Myr
  while t<tend-dt/2:
    print t        
    rad.evolve_model(t+dt)
    
    t+=dt
    i+=1

    if i%10==0: 
      p=sph.particles.copy()
      channel = rad.particles.new_channel_to(p)
      channel.copy_attributes(["xion","rho","u"])  
      print (mu/constants.kB*p.u.amin()).in_(units.K), (mu/constants.kB*p.u.amax()).in_(units.K)
      if write_snapshots:
        write_set_to_file(p,label+"-%6.6i"%i,"amuse",
                          append_to_file=False)    


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


def test(evolver, rad_parameters=None):    
  N=100000
  Ns=10
  L=15. | units.kpc
  dt=1. | units.Myr  
  tend=200.| units.Myr

  if rad_parameters==None:
    rad_parameters=dict(box_size=2*L)
  else:
    rad_parameters["box_size"]=2*L

  sph,rad=iliev_test_5(N,Ns,L,dt/2.,rad_parameters=rad_parameters)

  evolver(sph,rad,tend,dt)

  p=rad.particles.copy()
  channel = sph.particles.new_channel_to(p)
  channel.copy_attributes(["vx","vy","vz"])  

  r=(((p.x)**2+(p.y)**2+(p.z)**2) ** 0.5).value_in(units.kpc)
  v=((p.vx**2+p.vy**2+p.vz**2)**0.5).value_in(units.kms)
  cs=(p.u**0.5).value_in(units.kms)
  xion=p.xion.number
  t=mu/(1.+xion)/constants.kB*p.u
  t=t.value_in(units.K)
  dens=(p.rho).value_in(units.amu/units.cm**3)
  pres=(p.u*p.rho).value_in( 100*units.g/units.cm/units.s**2)
  mach=v/cs+1.e-10

  pyplot.figure( figsize=(8,6) ) 
  pyplot.plot( r,xion,'r.')
  pyplot.savefig('rx.png')

  pyplot.figure( figsize=(8,6) ) 
  pyplot.plot( r,t,'r.')
  pyplot.savefig('rt.png')

  i=int(tend.value_in(units.Myr))
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



if __name__=="__main__":
  global label
  import time
  
  
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
  test(radhydro_evolve)
  t2=time.time()
  print "test4: monochromatic, dynamic, no thermal solver",
  print t2-t1,"sec"

  t1=time.time()
  label="test5"
  test(radhydro_evolve_fake, dict(number_of_freq_bins=5, blackbody_spectrum_flag=1) )
  t2=time.time()
  print "test5: 5 freq, dynamic, no thermal solver",
  print t2-t1,"sec"
  
  
  t1=time.time()
  label="test6"
  test(radhydro_evolve, dict(number_of_freq_bins=5, blackbody_spectrum_flag=1,thermal_evolution_flag=1) )
  t2=time.time()
  print "test6: 5 freq, dynamic, thermal solver",
  print t2-t1,"sec"

