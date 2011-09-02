import numpy
from amuse.units import units
from amuse.units import nbody_system

from amuse.units import constants
from amuse.community.simplex.interface import SimpleX
from amuse.community.fi.interface import Fi

from amuse.support.data import core

from amuse.ext.evrard_test import uniform_random_unit_cube,uniform_unit_sphere

from matplotlib import pyplot

from amuse.io import write_set_to_file

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
  sph.parameters.pboxsize=2*L
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

  x=sph.particles.x.value_in(nbody_system.length)
  y=sph.particles.y.value_in(nbody_system.length)
  z=sph.particles.z.value_in(nbody_system.length)

  del sph  
  return x,y,z


def iliev_test_5_ic( N=10000,
                  Ns=10,
                  L=15. | units.kpc ):

  mp=rhoinit*(2*L)**3/N
 
#  x,y,z=uniform_random_unit_cube(N).make_xyz()
  x,y,z=glass(N,target_rms=0.025)
  
  p=core.Particles(N)
  p.x=L*x
  p.y=L*y
  p.z=L*z
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
                  dt=None):
  
  gas,sources=iliev_test_5_ic(N,Ns,L)                


  conv=nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
     
  sph=Fi(conv,use_gl=False,mode='periodic',redirection='none')   
  sph.initialize_code()

  sph.parameters.use_hydro_flag=True
  sph.parameters.radiation_flag=False
  sph.parameters.self_gravity_flag=False
  sph.parameters.gamma=1
  sph.parameters.isothermal_flag=True
  sph.parameters.integrate_entropy_flag=False
  sph.parameters.timestep=dt  
  sph.parameters.verbosity=0
  sph.parameters.pboxsize=2*L
  sph.commit_parameters()
  sph.gas_particles.add_particles(gas)
  sph.commit_particles()

#  sph.start_viewer()
         
  rad=SimpleX(number_of_workers = 1,redirection='none')
  rad.initialize_code()

  rad.parameters.box_size=2*L
  rad.parameters.hilbert_order=0

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
  p=pa.copy()
  p.u=fake_internal_energy_from_xion(p.xion)
  channel = p.new_channel_to(sys.particles)
  channel.copy_attributes(["u"])
  
def radhydro_evolve(sph,rad,tend,dt):

  i=0
  write_set_to_file(sph.gas_particles,"dump-%6.6i"%i,"amuse",
                        append_to_file=False)    

  t=0. | units.Myr
  while t<tend-dt/2:
    print t        
    print "sph1"
    sph.evolve_model(t+dt/2)    
    print "rad"
    update_pos_rho(rad,sph.gas_particles)
    rad.evolve_model(t+dt)
    update_u(sph,rad.particles)
    print "sph2"
    sph.evolve_model(t+dt)    
    t+=dt
    i+=1
    write_set_to_file(sph.gas_particles,"dump-%6.6i"%i,"amuse",
                        append_to_file=False)    
    
    


if __name__=="__main__":
  
  N=100000
  Ns=20
  L=15. | units.kpc
  dt=1. | units.Myr  
  tend=500.| units.Myr

  sph,rad=iliev_test_5(N,Ns,L,dt/2.)

  radhydro_evolve(sph,rad,tend,dt)


  p=sph.gas_particles.copy()
  channel = rad.particles.new_channel_to(p)
  channel.copy_attributes(["xion"])

  r=(((p.x)**2+(p.y)**2+(p.z)**2) ** 0.5).value_in(units.kpc)
  v=(((p.vx)**2+(p.vy)**2+(p.vz)**2) ** 0.5).value_in(units.km/units.s)
  xion=p.xion.number
  t=mu/constants.kB*fake_internal_energy_from_xion(xion)

  t=t.value_in(units.K)

  pyplot.figure(figsize=(8,6))
  pyplot.semilogy(xion,t,'r .')
  pyplot.savefig("test.png")

  pyplot.figure(figsize=(8,6))
  pyplot.plot(r,v,'r .')
  pyplot.savefig("testv.png")

