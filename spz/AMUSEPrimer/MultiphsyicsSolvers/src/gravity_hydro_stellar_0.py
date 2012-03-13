import math
import numpy
from amuse.lab import *
from amuse.ext.protodisk import ProtoPlanetaryDisk
from optparse import OptionParser
from amuse.couple import bridge
from amuse.ext.evrard_test import uniform_unit_sphere

SPH_code = Fi
from amuse import datamodel
feedback_efficiency = 0.0
total_feedback_energy = 0 | units.erg

class Boxed_sph(SPH_code):
  def __init__(self, *args, **kargs):
    SPH_code.__init__(self, *args, **kargs)
    self.escapers=datamodel.Particles(0)
  
  def evolve_model(self, *args, **kargs):
    print "Evolve hydro model"
    self.stopping_conditions.out_of_box_detection.enable()
    #outofbox=0.9*self.parameters.pboxsize/2
    outofbox=0.9*self.parameters.periodic_box_size/2
#    outofbox=0.9*0.5*2000. | units.AU
    print "A"
    self.parameters.stopping_conditions_out_of_box_size = outofbox
    print "B", outofbox
    self.overridden().evolve_model(*args,**kargs)
    escapers = Particles()
    print "find escapers:"
    while self.stopping_conditions.out_of_box_detection.is_set():
      local_escapers=self.particles.select_array(
        lambda x,y,z: (x**2+y**2+z**2 > outofbox**2), ["x","y","z"])
      escapers.add_particles(local_escapers)
    print "Number of escapers: ", len(escapers)
    if len(escapers)>0:
      self.escapers.add_particles(escapers)
      self.particles.remove_particles(escapers)
    self.overridden().evolve_model(*args,**kargs)

def lmech(star):
  lm=0.5*mdot_leitherer1992v_reimers(star)*v_terminal_teff(star)**2
  if(len(star) == 1 ):
    return lm[0]
  else:
    return lm  

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

class SSEplus(SSE):
  def __init__(self,**options):
    self.model_time=zero
    SSE.__init__(self,**options)

  def evolve_model(self,tend):
    print "Evolve SSE model: ", tend
    if not hasattr(self.particles,'Emech'):
      self.particles.Lmech=lmech(self.particles)
      self.particles.Emech=(0.| units.Myr)*self.particles.Lmech

    stellar_type=self.particles.stellar_type.copy()
    prev_lm=self.particles.Lmech.copy()

    ret=SSE.evolve_model(self,tend)
    print "SSE evolved.."
    if tend>self.model_time:
      print "supernova.."
      dt=tend-self.model_time    
      self.model_time=tend
      lm=lmech(self.particles)
      self.particles.Lmech=lm.copy()
      self.particles.Emech_last_feedback=self.particles.Emech
      self.particles.Emech=self.particles.Emech+dt*(prev_lm+lm)/2. + \
                          e_supernova(self.particles.stellar_type,stellar_type)    
      print "return star"
    return ret        

def main(Mstar=1, Ndisk=100, Mdisk= 0.001, Rmin=1, Rmax=100, t_end=10, n_steps=10, filename="nbody.hdf5"):
    Mstar = Mstar | units.MSun
    Mdisk = Mdisk | units.MSun
    Rmin = Rmin | units.AU
    Rmax = Rmax | units.AU
    t_end = t_end | units.yr

    mgas=Mdisk/Ndisk
    converter=nbody_system.nbody_to_si(Mstar, Rmin)
    disk = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter, 
                              densitypower=1.5, 
                              Rmin=Rmin.value_in(units.AU), 
                              Rmax=Rmax.value_in(units.AU),q_out=1.0,
                              discfraction=1.0).result
    disk.h_smooth= Rmin/Ndisk

    star=Particles(1)
    star[0].mass = Mstar
    star.position = [0, 0, 0] | units.AU
    star.velocity = [0, 0, 0] | units.kms
    #Adding Emech_last_feedback as new parameter to star
    star[0].Emech_last_feedback=(0.| units.erg)
    star[0].grav_mass=(0.| units.MSun)

    stellar = SSEplus()
    stellar.parameters.metallicity = 0.02
    stellar.particles.add_particle(star)
    stellar.commit_particles()
    time = 6.65 | units.Myr
#    time = 0.5*t_end
    stellar.evolve_model(time)
#    print stellar.particles
    channel_from_stellar_to_framework = stellar.particles.new_channel_to(star)
    channel_from_stellar_to_framework.copy_attributes(["mass", "Emech", "Lmech"])
    print "El=", stellar.particles.Emech, stellar.particles.Lmech

    hydro=Boxed_sph(converter)
    hydro.parameters.periodic_box_size = 100*Rmax
    hydro.parameters.use_hydro_flag=True
    hydro.parameters.radiation_flag=False
    hydro.parameters.self_gravity_flag=True
    hydro.parameters.gamma=1.
    hydro.parameters.isothermal_flag=True
    hydro.parameters.integrate_entropy_flag=False
    hydro.parameters.timestep=0.25 * t_end
    hydro.gas_particles.add_particles(disk)
    outofbox=hydro.parameters.periodic_box_size
    escapers=hydro.particles.select_array(
        lambda x,y,z: (x**2+y**2+z**2 > outofbox**2), ["x","y","z"])
    print " escapers out of box",len(escapers),outofbox.in_(units.parsec)
    hydro.particles.remove_particles(escapers)
    print "Nhydro=", len(hydro.particles)

    gravity = Hermite(converter)
    gravity.particles.add_particles(star)

    channel_from_gravity_to_framework = gravity.particles.new_channel_to(star)
    channel_from_hydro_to_framework = hydro.particles.new_channel_to(disk)
    channel_from_framework_to_gravity = star.new_channel_to(gravity.particles)

    moving_bodies = ParticlesSuperset([star, disk])
    moving_bodies.move_to_center()

    write_set_to_file(moving_bodies, filename, 'hdf5')
    
    gravity_hydro = bridge.Bridge(use_threading=False)
    gravity_hydro.add_system(gravity, (hydro,) )
    gravity_hydro.add_system(hydro, (gravity,) )

    Etot_init = gravity_hydro.kinetic_energy + gravity_hydro.potential_energy
    Etot_prev = Etot_init

    time = 0.0 | t_end.unit
    dt = t_end/float(n_steps)
    gravity_hydro.timestep = 0.25*dt
    while time < t_end:
        time += dt

        stellar.evolve_model(time)
        channel_from_stellar_to_framework.copy_attributes(["mass", "grav_mass"])
        channel_from_stellar_to_framework.copy_attributes(["Emech", "Lmech", "Emech_last_feedback"])
        print "El=", stellar.particles.mass, stellar.particles.grav_mass, stellar.particles.Emech, stellar.particles.Lmech, stellar.particles.Emech_last_feedback
        channel_from_framework_to_gravity.copy_attributes(["mass"])

#        print "sph:", disk
        gravity_hydro.evolve_model(time)

        Etot_prev_se = gravity_hydro.kinetic_energy + gravity_hydro.potential_energy

        channel_from_gravity_to_framework.copy()
        channel_from_hydro_to_framework.copy()
        write_set_to_file(moving_bodies, filename, 'hdf5')

#        mechanical_feedback(moving_bodies)
        channel_from_gravity_to_framework.copy_attributes(["x","y","z","vx","vy","vz"])
        new_sph=datamodel.Particles(0)
        star_particles = stellar.particles
        star_particles.dmass=star_particles.grav_mass-star_particles.mass
        star_particles.u=(star_particles.Emech-star_particles.Emech_last_feedback)/star_particles.dmass    
        if numpy.any((star_particles.Emech-star_particles.Emech_last_feedback) < zero):
          print "feedback error"
          raise Exception
          
        print "dmass=", star_particles.dmass, mgas.in_(units.MSun)
        star_particles[0].dmass = 2*mgas
        losers=star_particles.select_array(lambda x: x > mgas, ["dmass"])
        print "eject particles:", losers.position.value_in(units.parsec)
        while len(losers) > 0:            
          add=datamodel.Particles(len(losers))
          add.mass=mgas
          add.h_smooth=0. | units.parsec
          dx,dy,dz=uniform_unit_sphere(len(losers)).make_xyz()
          print losers.radius
          feedback_radius = 10.*losers.radius
          add.x=losers.x+feedback_radius*dx      
          add.y=losers.y+feedback_radius*dy
          add.z=losers.z+feedback_radius*dz
          add.vx=losers.vx      
          add.vy=losers.vy
          add.vz=losers.vz
          add.u=feedback_efficiency*losers.u
          losers.grav_mass-=mgas
          losers.Emech_last_feedback+=mgas*losers.u
          new_sph.add_particles(add)  
          losers=star_particles.select_array(lambda x,y: x-y > mgas, ["grav_mass","mass"])

        print "gas particles added:", len(new_sph)
        if len(new_sph) == 0:
          return      
        hydro.gas_particles.add_particles(new_sph)
#        total_feedback_energy=total_feedback_energy+(new_sph.mass*new_sph.u).sum()
        print "total_feedback_energy=", total_feedback_energy+(new_sph.mass*new_sph.u).sum()
        print "new_sph.u=", new_sph.u**0.5

        channel_from_stellar_to_framework.copy_attribute("grav_mass","mass")  

#        channel=star_particles.new_channel_to(star_particles_addition)
#        channel_from_framework_to_stellar.copy_attribute("Emech_last_feedback")

        Ekin = gravity_hydro.kinetic_energy 
        Epot = gravity_hydro.potential_energy
        Etot = Ekin + Epot
        print "T=", time, 
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        print "star=", star
        Etot_prev = Etot

    gravity_hydro.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="float", default = 10,
                      help="number of diagnostics time steps [10]")
    result.add_option("-f", dest="filename", default = "gravhydro.hdf5",
                      help="output filename [gravhydro.hdf5]")
    result.add_option("--Ndisk", dest="Ndisk", type="int",default = 100,
                      help="number of disk particles [100]")
    result.add_option("--Mstar", dest="Mstar", type="float",default = 30,
                      help="stellar mass [30] MStar")
    result.add_option("--Mdisk", dest="Mdisk", type="float",default = 0.001,
                      help="disk mass [0.001] MStar")
    result.add_option("--Rmin", dest="Rmin", type="float",default = 1.0,
                      help="minimal disk radius [1] in AU")
    result.add_option("--Rmax", dest="Rmax", type="float",default = 100,
                      help="maximal disk radius [150] in AU")
    result.add_option("-t", dest="t_end", type="float", default = 1,
                      help="end time of the simulation [1] year")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

