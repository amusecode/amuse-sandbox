import math
import numpy
from amuse.lab import *
from amuse.ext.protodisk import ProtoPlanetaryDisk
from optparse import OptionParser
from amuse.couple import bridge
from amuse import datamodel

def Roche_radius(a, M, m) :
    q = M/m;
    qcrt = q**(1./3.)
    qcrt2 = qcrt*qcrt
    Rl = a*0.49*qcrt2/(0.6*qcrt2 + math.log(1 + qcrt))
    return Rl

class Boxed_sph(Fi):
  def __init__(self, *args, **kargs):
    Fi.__init__(self, *args, **kargs)
    self.escapers=datamodel.Particles(0)
  
  def evolve_model(self, *args, **kargs):
    self.stopping_conditions.out_of_box_detection.enable()
    #outofbox=0.9*self.parameters.pboxsize/2
    outofbox=1. | units.AU
    self.parameters.stopping_conditions_out_of_box_size = outofbox
    self.overridden().evolve_model(*args,**kargs)
    while self.stopping_conditions.out_of_box_detection.is_set():
      escapers=self.particles.select_array(
        lambda x,y,z: (x**2+y**2+z**2 > outofbox**2), ["x","y","z"])
      print "***", len(escapers)
      if len(escapers)>0:
        self.escapers.add_particles(escapers)
        self.particles.remove_particles(escapers)
      self.overridden().evolve_model(*args,**kargs)

def main(Mstar=1, Ndisk=100, Mdisk= 0.01, Rmin=1, Rmax=100, t_end=10, n_steps=10, filename="nbody.hdf5"):
#    numpy.random.seed(111)
    Mstar = Mstar | units.MSun
    Mdisk = Mdisk | units.MSun
    Rmin = Rmin | units.RSun
    Rmax = Rmax | units.RSun
    t_end = 0.1*t_end | units.hour
    dt = t_end/float(n_steps)

#    White dwarf
    stars=Particles(2)
    stars[0].mass=0.8 | units.MSun
    stars[1].mass=0.4 | units.MSun
    Mstars = stars.mass.sum()
    print "Mstars=", stars.mass, Mstars

    stars[0].radius=0 | units.RSun
    stars[1].radius=0 | units.RSun

    stars[0].position = [0, 0, 0] | units.RSun
    stars[0].velocity = [0, 0, 0] | units.kms

    r = a = 1.35363 | units.RSun
    v_circ = (constants.G*stars.mass.sum()*(2./r - 1./a)).sqrt().value_in(units.kms)

#    RLOF donor
    stars[1].position = [a.value_in(units.RSun), 0, 0] | units.RSun
    stars[1].velocity = [0, v_circ, 0] | units.kms

    #Disk around WD
    Rmax = 0.1*0.43*a 
    Rmin = 0.1*Rmax
    Roche_wd = Roche_radius(a, stars[0].mass, stars[1].mass)
    print "Roche=", Roche_wd, Rmax
    Rmax = 0.5*Roche_wd/Rmin
    print "Roche=", Roche_wd, Rmax
#    converter=nbody_system.nbody_to_si(stars.mass.sum(), Rmin)
    converter=nbody_system.nbody_to_si(stars[0].mass, Rmin)
    disk = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter, densitypower=1.5, 
                              Rmin=10, Rmax=100.0, q_out=1.0,
                              discfraction=0.0001).result
    print "Diskmass=", disk.mass.sum().in_(units.MSun)
    disk.move_to_center()
    moving_bodies = ParticlesSuperset([stars, disk])
    moving_bodies.move_to_center()

#    disk = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter, densitypower=1.5,# 
#                              Rmin=Rmin.value_in(units.RSun), 
#                              Rmax=Rmax.value_in(units.RSun),q_out=1.0,
#                              discfraction=1.0).result
    disk.h_smooth= 0.1*Rmin/Ndisk

#    hydro=Fi(converter)
#    hydro=Boxed_sph(convert_nbody=converter,use_gl=False)

    hydro=Boxed_sph(channel_type="ibis", hostname="ij",convert_nbody=converter,use_gl=False)
#    hydror = Fi(channel_type="ibis", hostname="bergvliet")
 
    hydro.parameters.use_hydro_flag=True
    hydro.parameters.radiation_flag=False
    hydro.parameters.self_gravity_flag=True
    hydro.parameters.gamma=1.
    hydro.parameters.isothermal_flag=True
    hydro.parameters.integrate_entropy_flag=False
    hydro.parameters.timestep=dt/80.
    hydro.parameters.periodic_box_size=2 | units.AU

    hydro.gas_particles.add_particles(disk)

    gravity = Hermite(converter)
    gravity.particles.add_particles(stars)
#    gravity.particles.move_to_center()

    channel_from_gravity_to_framework = gravity.particles.new_channel_to(stars)
    channel_from_hydro_to_framework = hydro.particles.new_channel_to(disk)

    write_set_to_file(moving_bodies, filename, 'hdf5')
    
    gravity_hydro = bridge.Bridge(use_threading=False)
    gravity_hydro.add_system(gravity, (hydro,) )
    gravity_hydro.add_system(hydro, (gravity,) )

    Etot_init = gravity_hydro.kinetic_energy + gravity_hydro.potential_energy
    Etot_prev = Etot_init

    time = 0.0 | t_end.unit
#    dt = t_end/float(n_steps)
    gravity_hydro.timestep = dt/10.
    while time < t_end:
        time += dt
        gravity_hydro.evolve_model(time)

        Etot_prev_se = gravity_hydro.kinetic_energy + gravity_hydro.potential_energy

        channel_from_gravity_to_framework.copy()
        channel_from_hydro_to_framework.copy()
        write_set_to_file(moving_bodies, filename, 'hdf5')

        Ekin = gravity_hydro.kinetic_energy 
        Epot = gravity_hydro.potential_energy
        Etot = Ekin + Epot
        print "T=", time, 
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        print "star=", stars
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
    result.add_option("--Mstar", dest="Mstar", type="float",default = 1,
                      help="stellar mass [1] MStar")
    result.add_option("--Mdisk", dest="Mdisk", type="float",default = 0.001,
                      help="disk mass [0.001] MStar")
    result.add_option("--Rmin", dest="Rmin", type="float",default = 1.0,
                      help="minimal disk radius [1] in AU")
    result.add_option("--Rmax", dest="Rmax", type="float",default = 100,
                      help="maximal disk radius [150] in AU")
    result.add_option("-t", dest="t_end", type="float", default = 1.0,
                      help="end time of the simulation [1] year")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

