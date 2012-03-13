"""
   Evolve a single star in a molecular cloud.
   Allow accretion from the cloud onto the star.
"""
import sys
import numpy
from amuse.lab import *
from optparse import OptionParser
import amuse.datamodel as core
from amuse.couple import bridge
#from amuse.community.seba.interface import SeBa

def e_supernova(stellar_type,prev_stellar_type):
  i=numpy.where( (stellar_type>=13 | units.stellar_type) &
                 (stellar_type<=15 | units.stellar_type) &
                 ((prev_stellar_type<13 | units.stellar_type) |
                  (prev_stellar_type>15 | units.stellar_type)) )[0] 
  n=len(stellar_type)
  e=numpy.array([0.]*n)
  e[i]=1.e51
  if(n == 1 ):
    return (units.erg).new_quantity(e[0])
  else:
    return (units.erg).new_quantity(e)

    
class SSEplus(SSE):
  def __init__(self,**options):
    self.model_time=zero
    SSE.__init__(self,**options)

  def evolve_model(self,tend):
    if not hasattr(self.particles,'Emech'):
      self.particles.Lmech=lmech(self.particles)
      self.particles.Emech=(0.| units.Myr)*self.particles.Lmech

    stellar_type=self.particles.stellar_type.copy()
    prev_lm=self.particles.Lmech.copy()
    
    ret=SSE.evolve_model(self,tend)
    if tend>self.model_time:
      dt=tend-self.model_time    
      self.model_time=tend
      lm=lmech(self.particles)
      self.particles.Lmech=lm.copy()
      self.particles.Emech=self.particles.Emech+dt*(prev_lm+lm)/2. + \
                          e_supernova(self.particles.stellar_type,stellar_type)    
    return ret        

class BoxedFi(Fi):
  def __init__(self, *args, **kargs):
    Fi.__init__(self, *args, **kargs)
    self.escapers=core.Particles(0)
  
  def evolve_model(self, *args, **kargs):
    self.stopping_conditions.out_of_box_detection.enable()
    outofbox=0.9*self.parameters.periodic_box_size/2
    self.parameters.stopping_conditions_out_of_box_size = outofbox
#    Fi.evolve_model(self,*args,**kargs)
    self.overridden().evolve_model(*args,**kargs)
    while self.stopping_conditions.out_of_box_detection.is_set():
      escapers=self.particles.select_array(
        lambda x,y,z: (x**2+y**2+z**2 > outofbox**2), ["x","y","z"])
      print "***", len(escapers)
      if len(escapers)>0:
        self.escapers.add_particles(escapers)
        self.particles.remove_particles(escapers)
#      Fi.evolve_model(self,*args, **kargs)
      self.overridden().evolve_model(*args,**kargs)

def create_binary(interface, stars,ecc=0,a=1| units.AU):

  m1 = stars[0].mass
  m2 = stars[1].mass
  r1 = stars[0].radius
  r2 = stars[1].radius
  mu=constants.G*(m1+m2)
  f1=m2/(m1+m2)
  f2=m1/(m1+m2)
  rmax=a*(1+ecc)
  r0=rmax
  print a.in_(units.AU),r0.in_(units.AU)
  
  h=(a*mu*(1-ecc**2))**0.5
  v0=h/r0

  bin=stars

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

  convert=nbody_system.nbody_to_si(bin.mass.sum(),a)
#  nb = interface(convert,redirection="none")
  nb = interface(convert)
  nb.particles.add_particles(bin)
  nb.particles.move_to_center()

  return nb

def mdot_leitherer1992_variant(star):
  m_msun=star.mass.value_in(units.MSun)
  l_lsun=star.luminosity.value_in(units.LSun)
  temp=star.temperature.value_in(units.K)
  z_zsun=1.
  return (units.MSun/units.yr).new_quantity(
    10**(-11.35+2.45*numpy.log10(l_lsun) - 1.1*numpy.log10(m_msun) +\
    -1.31*numpy.log10(temp)+0.80*numpy.log10(z_zsun)) )

def mdot_leitherer1992v_reimers(star):
  m_msun=star.mass.value_in(units.MSun)
  l_lsun=star.luminosity.value_in(units.LSun)
  r_rsun=star.radius.value_in(units.RSun)
  temp=star.temperature.value_in(units.K).clip(1000,8.e4)
  z_zsun=1.
  mdot=(10**(-12.+2.45*numpy.log10(l_lsun) - 1.1*numpy.log10(m_msun) +\
      -1.3*numpy.log10(temp)+0.80*numpy.log10(z_zsun)))
  i=numpy.where(temp<8000.)[0]
  mdot[i]=4.e-13*l_lsun[i]*r_rsun[i]/m_msun[i]
  return (units.MSun/units.yr).new_quantity(mdot)


def v_terminal_teff(star):
  t4=numpy.log10(star.temperature.value_in(units.K))-4.
  t4=t4.clip(0.,1.)
  return (30 | units.km/units.s) + ((4000 | units.km/units.s)*t4)

def e_supernova(star):
  i=numpy.where( (star.stellar_type>=13 | units.stellar_type) &
    (star.stellar_type<=15 | units.stellar_type) )[0] 
  n=len(star)
  e=numpy.array([0.]*n)
  e[i]=1.e51
  if(n == 1 ):
    return (units.erg).new_quantity(e[0])
  else:
    return (units.erg).new_quantity(e)


def lmech(star):
  lm=0.5*mdot_leitherer1992v_reimers(star)*v_terminal_teff(star)**2
  if(len(star) == 1 ):
    return lm[0]
  else:
    return lm  

def mechanical_feedback(stellar, gravity_hydro):
    
    star_particles=stellar.particles.copy_to_memory()
    
    channel= self.particles.new_channel_to(star_particles)

    channel.copy_attributes(["x","y","z","vx","vy","vz"])
    channel.copy_attribute("mass","grav_mass")
    del channel

    channel=self.star_particles_addition.new_channel_to(star_particles)
    channel.copy_attribute("Emech_last_feedback")
    del channel

    new_sph=datamodel.Particles(0)
    
    star_particles.dmass=star_particles.grav_mass-star_particles.mass
    star_particles.u=(star_particles.Emech-star_particles.Emech_last_feedback)/star_particles.dmass    
    if numpy.any((star_particles.Emech-star_particles.Emech_last_feedback) < zero):
      print "feedback error"
      raise Exception
          
    losers=star_particles.select_array(lambda x: x > self.mgas, ["dmass"])
    print losers.x.value_in(units.parsec)
    print losers.y.value_in(units.parsec)   
    print losers.z.value_in(units.parsec)   
    while len(losers) > 0:            
      add=datamodel.Particles(len(losers))
      add.mass=self.mgas
      add.h_smooth=0. | units.parsec
      dx,dy,dz=uniform_unit_sphere(len(losers)).make_xyz()
      add.x=losers.x+self.feedback_radius*dx      
      add.y=losers.y+self.feedback_radius*dy
      add.z=losers.z+self.feedback_radius*dz
      add.vx=losers.vx      
      add.vy=losers.vy
      add.vz=losers.vz
      add.u=self.feedback_efficiency*losers.u
      losers.grav_mass-=self.mgas
      losers.Emech_last_feedback+=self.mgas*losers.u
      new_sph.add_particles(add)  
      losers=star_particles.select_array(lambda x,y: x-y > self.mgas, ["grav_mass","mass"])


    print "gas particles added:", len(new_sph)
    if len(new_sph) == 0:
      return      
    self.sph.gas_particles.add_particles(new_sph)
    self.total_feedback_energy=self.total_feedback_energy+(new_sph.mass*new_sph.u).sum()
    print new_sph.u**0.5

    channel = star_particles.new_channel_to(self.particles)
    channel.copy_attribute("grav_mass","mass")  
    del channel

    channel=star_particles.new_channel_to(self.star_particles_addition)
    channel.copy_attribute("Emech_last_feedback")
    del channel

def main(Mstar=1.0, z=0.02, Mmc=100, Nsph=1000, Rmc=1.0, t_end=10, filename="stellar.hdf5"):
    t_end = t_end | units.Myr
    Mstar = Mstar | units.MSun
    Mmc = Mmc | units.MSun
    Rmc = Rmc  | units.parsec
    stellarfile = filename
    hydrofile = "hydro.hdf5"

    # set up the star
    stellar = SSEplus()
    stellar.parameters.metallicity = z
    m1=Mstar
    m2=Mstar
    stars = Particles(2)
    stars[0].mass = m1
    stars[1].mass = m2
    stellar.particles.add_particle(stars)
    stellar.commit_particles()
    stellar.evolve_model(1.0|units.yr)
    channel_from_stellar_to_framework = stellar.particles.new_channel_to(stars)
    print "stars=", stellar.particles

    gravity=create_binary(ph4,stellar.particles,ecc=0.15944,a=1.0 | units.AU)

    channel_from_gravity_to_framework = gravity.particles.new_channel_to(stars)
    channel_from_framework_to_gravity = stars.new_channel_to(gravity.particles)
    channel_from_stellar_to_framework.copy_attributes(["mass","radius"])

    # set up the hydro cloud
    converter=nbody_system.nbody_to_si(Mmc, Rmc)
    cloud = new_plummer_gas_model(Nsph, convert_nbody=converter)
    outofbox=10*Rmc
    escapers=cloud.select_array(
        lambda x,y,z: (x**2+y**2+z**2 > outofbox**2), ["x","y","z"])
    print "**",len(escapers),outofbox.in_(units.parsec)
    cloud.remove_particles(escapers)
    print len(cloud)

    # initialize hydro code
    hydro=BoxedFi(convert_nbody=converter,use_gl=False)
    hydro.parameters.periodic_box_size = 10*Rmc
    hydro.parameters.timestep=converter.to_si(1|nbody_system.time)/32
    hydro.parameters.self_gravity_flag=False
    hydro.gas_particles.add_particles(cloud)
    channel_from_hydro_to_framework = hydro.gas_particles.new_channel_to(cloud)

    write_set_to_file(gravity.particles, stellarfile, 'hdf5')
    write_set_to_file(hydro.particles, hydrofile, 'hdf5')

    gravity_hydro = bridge.Bridge(use_threading=False)
    gravity_hydro.add_system(gravity, (hydro,), False )
    gravity_hydro.add_system(hydro, (gravity,), False )
    gravity_hydro.timestep = t_end/10.

    Mtot_init = stellar.particles.mass.sum()
    time = 0.0 | t_end.unit
    dt = t_end/100
    while time < t_end:
        time += dt
        stellar.evolve_model(time)
#        time = stellar.model_time

        channel_from_stellar_to_framework.copy_attributes(["mass","radius"])
        channel_from_framework_to_gravity.copy_attributes(["mass"])

        gravity_hydro.evolve_model(time)
        channel_from_gravity_to_framework.copy()
        channel_from_hydro_to_framework.copy()
 
        write_set_to_file(stellar.particles, stellarfile, 'hdf5')
        write_set_to_file(hydro.particles, hydrofile, 'hdf5')

#        mechanical_feedback(stellar, gravity_hydro)

        print "before sink:", stellar.particles
        print "Ngas=", len(hydro.particles)
#        sink_particles(gravity.particles,hydro.particles,Raccretion=100000.0)
        print "after sink:", stellar.particles
        print "Ngas=", len(hydro.particles)

    gravity.stop()
    hydror.stop()
    stellar.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "stellar.hdf5",
                      help="output filename [stellar.hdf5]")
    result.add_option("--Mstar", dest="Mstar", type="float",default = 1,
                      help="stellar mass [1] MSun")
    result.add_option("--Mmc", dest="Mmc", type="float",default = 100,
                      help="molecular cloud mass [100] MSun")
    result.add_option("--Rmc", dest="Rmc", type="float",default = 100,
                      help="molecular cloud radius [100] pc")
    result.add_option("--Nsph", dest="Nsph", type="int",default = 1000,
                      help="number of sph particles in cloud [1000]")
    result.add_option("-t", dest="t_end", type="float", default = 1.0,
                      help="end time of the simulation [1] Myr")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

