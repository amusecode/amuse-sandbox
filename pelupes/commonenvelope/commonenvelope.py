import numpy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot


from amuse.units import units,constants
from amuse.community.sse.interface import SSE
from amuse.community.evtwin.interface import EVtwin
from amuse.community.mesa.interface import MESA
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH,pickle_stellar_model
from amuse.plot import sph_particles_plot, native_plot
from amuse.datamodel import Particle,Particles

from amuse.community.fi.interface import Fi
from amuse.community.gadget2.interface import Gadget2

from amuse.units import nbody_system

#import logging
#logging.basicConfig(level=logging.DEBUG)

# returns r1/a, q=m1/m2
def roche_radius_over_a(q):
  return 0.49* q**(2./3.)/ (0.6*q**(2./3.) + numpy.log(1.+q**(1./3.)))



def generate_stellar_model(mass=1.05|units.MSun,filename='./star.pkl',Rtarget=2.2e12 | units.cm,time=None):
  
    star =  Particle()
    star.mass = mass
    
    stellar_evolution = MESA()
    print stellar_evolution.parameters
    se_star = stellar_evolution.particles.add_particle(star)
    
    print "Evolving", star.mass, "star with", stellar_evolution.__class__.__name__, "up to", 100 | units.Myr

    tend=500 | units.Myr
    dt=10. | units.Myr
    t=0.*tend
    while se_star.radius<Rtarget:
      stellar_evolution.evolve_model()
      t=stellar_evolution.model_time
      print t.in_(units.Myr),se_star.radius

    print se_star  

    pickle_stellar_model(se_star,filename)

    stellar_evolution.stop()

def sse_evolve(mass=1.05| units.MSun,tend=9990. | units.Myr,Rtarget=2.2e12 | units.cm):
    star =  Particle()
    star.mass = mass
    
    stellar_evolution = SSE()
    se_star = stellar_evolution.particles.add_particle(star)

    while se_star.radius<Rtarget:
      stellar_evolution.evolve_model()
      t=stellar_evolution.model_time
      print t.in_(units.Myr),se_star.radius

    print se_star  


def create_particles(N=10000,core_mass=0.32 | units.MSun,filename="./star.pkl"):
    model = convert_stellar_model_to_SPH(
        None, 
        N,
        with_core_particle=True,
        target_core_mass=core_mass,
        pickle_file=filename#,base_grid_options=dict(type="random")
    )
    return model.core_particle.as_set(),model.gas_particles

def make_map(sph,L,N=200):

    x,y=numpy.indices( ( N+1,N+1 ))

    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=x*0.
        
    vx=units.kms(numpy.zeros_like(x.number))
    vy=units.kms(numpy.zeros_like(x.number))
    vz=units.kms(numpy.zeros_like(x.number))

    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
#    kin_energy=0.5*(rhovx**2+rhovy**2+rhovz**2)/rho
#    internal_energy=abs(rhoe-kin_energy)/rho
#    temp=global_mu / constants.kB * internal_energy

    rho=rho.reshape((N+1,N+1))
#    temp=temp.reshape((N+1,N+1))

    return numpy.transpose(rho)#,numpy.transpose(temp)

def output_map(sph,i,L=60| units.RSun,N=200):
      rho=make_map(sph,L,N=N)
      
      f=pyplot.figure(figsize=(8,8))
      LL=L.number/2

      vmax=(sph.gas_particles.rho.amax()).value_in(units.amu/units.cm**3)
      vmin=(sph.gas_particles.rho.amin()).value_in(units.amu/units.cm**3)
      vmax=2.5e20
      vmin=7.e17
      
      pyplot.imshow(numpy.log10(rho.value_in(units.amu/units.cm**3)+vmin/100),
          extent=[-LL,LL,-LL,LL],vmin=numpy.log10(vmin),vmax=numpy.log10(vmax),origin='lower')
      pyplot.xlabel(L.unit)
      pyplot.savefig('map-%6.6i.png'%i)
      f.clear()
      pyplot.close(f)

def output_map(sph,i,L=60| units.RSun):
      
      f=pyplot.figure(figsize=(8,8))
      LL=L.number/2

      xcm,ycm,zcm=sph.particles.center_of_mass().value_in(L.unit)

      x=sph.gas_particles.x.value_in(L.unit)-xcm
      y=sph.gas_particles.y.value_in(L.unit)-ycm
      
      pyplot.plot(x,y,'r.')

      x=sph.dm_particles.x.value_in(L.unit)-xcm
      y=sph.dm_particles.y.value_in(L.unit)-ycm
      pyplot.plot(x,y,'b+',markersize=12.,mew=2.)

      pyplot.xlim(-LL,LL)
      pyplot.ylim(-LL,LL)
      pyplot.savefig('pmap-%6.6i.png'%i)

      f.clear()
      pyplot.close(f)

def radial_u_plot(sph,i,L=60| units.RSun):
      
      f=pyplot.figure(figsize=(8,8))
      LL=L.number/2

      cx=sph.dm_particles[0].x
      cy=sph.dm_particles[0].y
      cz=sph.dm_particles[0].z

      r=((sph.gas_particles.x-cx)**2+(sph.gas_particles.y-cy)**2+(sph.gas_particles.z-cz)**2)**0.5
      r=r.value_in(units.RSun)
      u=sph.gas_particles.u.value_in(units.kms**2)

      pyplot.plot(r,u,'r.')


      pyplot.xlim(0,LL*2)
      pyplot.ylim(0,3.e4)
      pyplot.savefig('umap-%6.6i.png'%i)

      f.clear()
      pyplot.close(f)

def radial_rho_plot(sph,i,L=60| units.RSun):
      
      f=pyplot.figure(figsize=(8,8))
      LL=L.number/2

      cx=sph.dm_particles[0].x
      cy=sph.dm_particles[0].y
      cz=sph.dm_particles[0].z

      r=((sph.gas_particles.x-cx)**2+(sph.gas_particles.y-cy)**2+(sph.gas_particles.z-cz)**2)**0.5
      r=r.value_in(units.RSun)
      rho=sph.gas_particles.rho.value_in(units.amu/units.cm**3)

      pyplot.plot(r,rho,'r.')


      pyplot.xlim(0,LL*2)
      pyplot.ylim(0,1.e21)
      pyplot.savefig('rhomap-%6.6i.png'%i)

      f.clear()
      pyplot.close(f)


def radial_A_plot(sph,i,L=60| units.RSun):
      
      f=pyplot.figure(figsize=(8,8))
      LL=L.number/2

      cx=sph.dm_particles[0].x
      cy=sph.dm_particles[0].y
      cz=sph.dm_particles[0].z

      r=((sph.gas_particles.x-cx)**2+(sph.gas_particles.y-cy)**2+(sph.gas_particles.z-cz)**2)**0.5
      r=r.value_in(units.RSun)
      u=sph.gas_particles.u.value_in(units.kms**2)
      rho=sph.gas_particles.rho.value_in(units.amu/units.cm**3)

      gamma=5./3.
      A=(gamma-1)*rho*u/rho**gamma

      pyplot.plot(r,A,'r.')

      pyplot.xlim(0,LL*2)
      pyplot.ylim(0,6.e-10)
      pyplot.savefig('Amap-%6.6i.png'%i)

      f.clear()
      pyplot.close(f)


def eos_plot(sph,i):
      
      f=pyplot.figure(figsize=(8,8))

      rho=sph.gas_particles.rho.value_in(units.amu/units.cm**3)
      u=sph.gas_particles.u.value_in(units.kms**2)
      p=rho*u

      j=len(p)/2
      pyplot.loglog(rho,p[j]*(rho/rho[j])**1.6667,'g.')
      pyplot.loglog(rho,p[j]*(rho/rho[j])**1.,'b.')
      pyplot.loglog(rho,p,'r.')

      pyplot.savefig('eos-%6.6i.png'%i)

      f.clear()
      pyplot.close(f)


def sph_star(N=1000,core_mass=0.32 | units.MSun,code=Fi):
    core,gas=create_particles(N=N,core_mass=core_mass)
            
    dt=min((1|units.RSun)/(gas.u)**0.5)/2
    
    print dt.in_(units.hour)

    UnitLength=1.| units.RSun
    UnitMass=.001 | units.MSun 
  
    convert=nbody_system.nbody_to_si(UnitLength,UnitMass)

    print convert.to_nbody(gas.u.amin())
    print convert.to_nbody(gas.u.amax())

    sph=code(convert,number_of_workers=2)#,redirection="none")#number_of_workers=1,use_gl=False,debugger='gdb')
  
    sph.parameters.periodic_box_size=1000. | units.RSun

    if code==Fi:
      sph.parameters.timestep=dt
  
    if code==Gadget2:
      sph.parameters.max_size_timestep=dt
      sph.parameters.time_max=240 | units.day
      sph.parameters.epsilon_squared=(2*core[0].radius)**2
      sph.parameters.n_smooth=64
    
#    print sph.parameters
    core.radius*=2
    print core.radius.in_(units.RSun)
    
    sph.gas_particles.add_particles(gas)
    sph.dm_particles.add_particles(core)

    print convert.to_nbody(sph.gas_particles.rho.amin())
    print convert.to_nbody(sph.gas_particles.rho.amax())

    return sph,dt

def write_model(sph,i,label=""):
    from amuse.io import write_set_to_file

    print sph.gas_particles[-1].u.in_(units.kms**2)
    
    dm_particles=sph.dm_particles.copy()
    gas_particles=sph.gas_particles.copy()

    print gas_particles[-1].u.in_(units.kms**2)
    
    cm=dm_particles.center_of_mass()
    cmv=dm_particles.center_of_mass_velocity()
    
    dm_particles.position-=cm
    dm_particles.velocity-=cmv
    gas_particles.position-=cm
    gas_particles.velocity-=cmv
        
    print gas_particles[-1].u.in_(units.kms**2)
        
    write_set_to_file(dm_particles,label+'core-%8.8i.amuse'%i,'amuse')
    write_set_to_file(gas_particles,label+'gas-%8.8i.amuse'%i,'amuse')

def single_star_test(N=1000,core_mass=0.32 | units.MSun,tend=60| units.day):
    sph,dtsph=sph_star(N=N,core_mass=core_mass)
    
    print (sph.gas_particles.rho.amax()).in_(units.amu/units.cm**3), 
    print (sph.gas_particles.rho.amin()).in_(units.amu/units.cm**3)

    print len(sph.gas_particles)

    output_map(sph,0)
    radial_A_plot(sph,0)
    radial_u_plot(sph,0)
    radial_rho_plot(sph,0)
#    eos_plot(sph,0)
    
    r=(sph.gas_particles.x**2+sph.gas_particles.y**2+sph.gas_particles.z**2)**0.5
    

    print "rmax:",r.amax().in_(units.RSun)
    v=(sph.gas_particles.vx**2+sph.gas_particles.vy**2+sph.gas_particles.vz**2)**0.5
    print "vmax:",v.amax().in_(units.kms)

    print "hmin:",sph.gas_particles.h_smooth.amin().value_in(units.RSun)

#    print sph.dm_particles.position/(1|units.RSun)
#    print sph.dm_particles.velocity/(1|units.kms)

#    print sph.dm_particles[0].radius.in_(units.RSun)

    dt=int((4| units.hour)/dtsph) *dtsph

    i=0
    t=0.*tend
    while t< tend-dt/2:
      t=t+dt
      i=i+1
      sph.evolve_model(t)

      r=(sph.gas_particles.x**2+sph.gas_particles.y**2+sph.gas_particles.z**2)**0.5
      v=(sph.gas_particles.vx**2+sph.gas_particles.vy**2+sph.gas_particles.vz**2)**0.5
      print "t,rmax,vmax:",t.in_(units.day),r.amax().in_(units.RSun),v.amax().in_(units.kms)

 #     print sph.dm_particles.position/(1|units.RSun)
 #     print sph.dm_particles.velocity/(1|units.kms)
      
      output_map(sph,i)
      radial_u_plot(sph,i)
      radial_A_plot(sph,i)
      radial_rho_plot(sph,i)
      
      r=(sph.gas_particles.x**2+sph.gas_particles.y**2+sph.gas_particles.z**2)**0.5
      a=numpy.where(r> 50 | units.RSun)[0]
      print "escape:",len(a)
      
      v=(sph.gas_particles.vx**2+sph.gas_particles.vy**2+sph.gas_particles.vz**2)**0.5
      a=numpy.where(v> 30 | units.kms)[0]
      print "highv:",len(a)
    write_model(sph,N)  

def rotating_star_test(N=1000,core_mass=0.32 | units.MSun,tend=60| units.day
                       ,Redge=2.2e12| units.cm,V_Vcrit=0.95,a=4.4e12| units.cm,m2=0.6 | units.MSun,
                       code=Fi):
    sph,dtsph=sph_star(N=N,core_mass=core_mass,code=code)

    m1=sph.particles.total_mass()
    
    print m1,m2

    mu=constants.G*(m1+m2)
    P=2*numpy.pi*a**(3./2)/mu**0.5

    gas=sph.gas_particles.copy()

    # uniform rotation
    Vrot=2*numpy.pi*Redge/P*V_Vcrit
    r=(gas.x**2+gas.y**2)**0.5
    gas.vy+=-Vrot*gas.x/Redge
    gas.vx+=Vrot*gas.y/Redge

    channel=gas.new_channel_to(sph.gas_particles)
    channel.copy_attributes(["vx","vy","vz"])
    
    output_map(sph,0)
    radial_A_plot(sph,0)
    radial_u_plot(sph,0)
    radial_rho_plot(sph,0)
#    eos_plot(sph,0)
    
    dt=int((4| units.hour)/dtsph) *dtsph

    i=0
    t=0.*tend
    while t< tend-dt/2:
      t=t+dt
      i=i+1
      sph.evolve_model(t)

      print "t:",t.in_(units.day)
      
      output_map(sph,i)
      radial_u_plot(sph,i)
      radial_A_plot(sph,i)
      radial_rho_plot(sph,i)
      
    write_model(sph,N)  

    
def setup_ce(code=Fi,N=1000,Redge=2.2e12| units.cm ,a=4.3e12 | units.cm,
             m2=0.6| units.MSun,ecc=0):
    from amuse.io import read_set_from_file
    
    core=read_set_from_file('core-%8.8i.amuse'%N,'amuse')
    gas=read_set_from_file('gas-%8.8i.amuse'%N,'amuse')

    secondary=Particles(1)    
    secondary.radius=core.radius
    secondary.mass=m2
        
    m1=core.total_mass()+gas.total_mass()
    m2=secondary.total_mass()

    print m1,m2

    mu=constants.G*(m1+m2)
    P=2*numpy.pi*a**(3./2)/mu**0.5
    print P.in_(units.day)

    f1=m2/(m1+m2)
    f2=m1/(m1+m2)

    rmax=a*(1+ecc)
    r0=rmax
    print a.in_(units.AU),r0.in_(units.AU)
  
    h=(a*mu/(1-ecc**2))**0.5
    v0=h/r0

    secondary.x=f2*r0
    secondary.y=0.*r0
    secondary.z=0.*r0
    secondary.vy=-f2*v0
    secondary.vx=0.*v0
    secondary.vz=0.*v0

    core.x+=-f1*r0
    core.vy+=f1*v0
    gas.x+=-f1*r0
    gas.vy+=f1*v0
 
    dt=min((1|units.RSun)/(gas.u)**0.5)/2
    
    UnitLength=1.| units.RSun
    UnitMass=.001 | units.MSun 
  
    convert=nbody_system.nbody_to_si(UnitLength,UnitMass)
    
    if code==Fi:
      sph=code(convert)
      sph.parameters.periodic_box_size=10000. | units.RSun
      sph.parameters.timestep=dt
  
    if code==Gadget2:
      sph=code(convert,number_of_workers=4)
      sph.parameters.periodic_box_size=10000. | units.RSun
      sph.parameters.max_size_timestep=dt
      sph.parameters.time_max=1000. | units.day
      sph.parameters.epsilon_squared=(core[0].radius)**2
      sph.parameters.n_smooth=64

    sph.gas_particles.add_particles(gas)
    sph.dm_particles.add_particles(core)
    sph.dm_particles.add_particles(secondary)

    return sph,dt
    
def ce_evolve(code=Fi,N=1000,Redge=2.2e12| units.cm,a=4.4e12| units.cm,m2=0.6 | units.MSun,
                outstep=100,tend=120.| units.day):
    sph,dtsph=setup_ce(code=code,N=N,Redge=Redge,a=a,m2=m2)
    
    output_map(sph,0,L=200| units.RSun)
    radial_A_plot(sph,0)
    radial_u_plot(sph,0)
    radial_rho_plot(sph,0)
    
    dt=int((4| units.hour)/dtsph) *dtsph

    i=0
    t=0.*tend
    while t< tend-dt/2:
      t=t+dt
      i=i+1
      sph.evolve_model(t)
      print "t:",t.in_(units.day)      

      output_map(sph,i,L=200| units.RSun)
      radial_u_plot(sph,i)
      radial_A_plot(sph,i)
      radial_rho_plot(sph,i)
      
      if i%outstep==0:
        write_model(sph,i,label="ce-")  

    write_model(sph,i,label="ce-")  

      
if __name__ in ("__main__", "__plot__"):
    m1=3.17| units.MSun
    m2=0.57| units.MSun
    r_a=roche_radius_over_a(m1/m2)
    a=71.8| units.RSun
    print r_a,(r_a*a).in_(units.RSun)
    mcore=(0.56 | units.MSun)
  
#    ce_evolve(code=Gadget2,N=10000,m2=m2,Redge=r_a*a,a=0.95*a,tend=180.| units.day)

    generate_stellar_model(mass=m1,Rtarget=r_a*a)
#    sse_evolve(mass=m1,Rtarget=r_a*a)

    single_star_test(N=1000,tend=120.|units.day,core_mass=mcore)
#    rotating_star_test(N=10000,tend=120.|units.day,core_mass=mcore,m2=m2,Redge=r_a*a,a=a,
#      code=Gadget2,V_Vcrit=0.95)

#    import cProfile
#    cProfile.run("sph_particles = create_particles()","prof_100k")
#    sph_particles = create_particles(N=100000)
#    print "done"
#    native_plot.figure(figsize = (6, 6), dpi = 100)
#    sph_particles_plot(sph_particles)
#    native_plot.savefig("test.png")

#print roche_radius_over_a(1.05/0.6)*4.3
