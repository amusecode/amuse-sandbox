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

def sink_particles(sinks,sources,FRaccretion=0.1 | units.AU):

#  R2=Raccretion**2
  all_lost_particles = Particles()
  for s in sinks:
     xs,ys,zs=s.x,s.y,s.z
     R2 = (FRaccretion*s.radius)**2
     insink=sources.select_array(lambda x,y,z: (x-xs)**2+(y-ys)**2+(z-zs)**2 < R2,['x','y','z'])  
     print "Racc=", R2.sqrt().as_quantity_in(units.AU), len(insink)
     if len(insink) > 0:
       cm=s.position*s.mass
       p=s.velocity*s.mass
       s.mass+=insink.total_mass()
       s.position=(cm+insink.center_of_mass()*insink.total_mass())/s.mass
       s.velocity=(p+insink.total_momentum())/s.mass
       print "Nlost=", len(insink), "of", len(sources)
#       print insink[0:1]   
       all_lost_particles.add_particles(insink)
  if len(all_lost_particles)>0:
    sources.remove_particles(all_lost_particles)
    print "left over N=", len(sources)

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

  convert=nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
#  nb = interface(convert,redirection="none")
  nb = interface(convert)
  nb.particles.add_particles(bin)
  return nb

def main(Mstar=1.0, z=0.02, Mmc=100, Nsph=1000, Rmc=1.0, t_end=10, filename="stellar.hdf5"):
    t_end = t_end | units.Myr
    Mstar = Mstar | units.MSun
    Mmc = Mmc | units.MSun
    Rmc = Rmc  | units.parsec
    stellarfile = filename
    hydrofile = "hydro.hdf5"

    # set up the star
    stellar = SSE()
    stellar.parameters.metallicity = z
    m1=Mstar
    m2=0.5*Mstar
    stars = Particles(2)
    stars[0].mass = m1
    stars[1].mass = m2
    stellar.particles.add_particle(stars)
    stellar.commit_particles()
    stellar.evolve_model(1.0|units.yr)
    channel_from_stellar_to_framework = stellar.particles.new_channel_to(stars)
    print "stars=", stellar.particles

    gravity=create_binary(Hermite,stellar.particles,ecc=0.15944,a=1.0 | units.AU)
    channel_from_gravity_to_framework = gravity.particles.new_channel_to(stars)
    channel_from_framework_to_gravity = stars.new_channel_to(gravity.particles)
    channel_from_stellar_to_framework.copy_attributes(["mass","radius"])

    # set up the hydro cloud
    converter=nbody_system.nbody_to_si(Mmc, Rmc)
    cloud = new_plummer_gas_model(Nsph, convert_nbody=converter)
    outofbox=10. | units.parsec
    escapers=cloud.select_array(
        lambda x,y,z: (x**2+y**2+z**2 > outofbox**2), ["x","y","z"])
    print "**",len(escapers),outofbox.in_(units.parsec)
    cloud.remove_particles(escapers)
    print len(cloud)

    # initialize hydro code
    hydro=BoxedFi(convert_nbody=converter,use_gl=False)
    hydro.parameters.periodic_box_size=20. | units.parsec
    hydro.parameters.timestep=converter.to_si(1|nbody_system.time)/32
    hydro.parameters.self_gravity_flag=False
    hydro.gas_particles.add_particles(cloud)
    hydro.commit_particles()
    channel_from_hydro_to_framework = hydro.gas_particles.new_channel_to(cloud)

    write_set_to_file(gravity.particles, stellarfile, 'hdf5')
    write_set_to_file(hydro.particles, hydrofile, 'hdf5')

    gravity_hydro = bridge.Bridge(use_threading=False)
    gravity_hydro.add_system(gravity, (hydro,), False )
    gravity_hydro.add_system(hydro, (gravity,), False )
    gravity_hydro.timestep = t_end/10.

    Mtot_init = stellar.particles.mass.sum()
    time = 0.0 | t_end.unit
    while time < t_end:

        print "before stellar:", stellar.particles
        stellar.evolve_model()
        print "after stellar:", stellar.particles
        time = stellar.model_time
        print "time=", time

        channel_from_stellar_to_framework.copy_attributes(["mass","radius"])
        channel_from_framework_to_gravity.copy_attributes(["mass"])

        gravity_hydro.evolve_model(time)
#        hydro.evolve_model(time)
        channel_from_gravity_to_framework.copy()
        channel_from_hydro_to_framework.copy()
 
        write_set_to_file(stellar.particles, stellarfile, 'hdf5')
        write_set_to_file(hydro.particles, hydrofile, 'hdf5')

        print "before sink:", stellar.particles
        print "Ngas=", len(hydro.particles)
        sink_particles(gravity.particles,hydro.particles,FRaccretion=1000000.0)
        print "after sink:", stellar.particles
        print "Ngas=", len(hydro.particles)

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
    result.add_option("--Rmc", dest="Rmc", type="float",default = 1,
                      help="molecular cloud radius [10] pc")
    result.add_option("--Nsph", dest="Nsph", type="int",default = 1000,
                      help="number of sph particles in cloud [1000]")
    result.add_option("-t", dest="t_end", type="float", default = 100.0,
                      help="end time of the simulation [100] Myr")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

