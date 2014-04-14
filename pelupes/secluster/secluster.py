import numpy 
import cPickle

from matplotlib import pyplot
from amuse.units import units
from amuse.units import nbody_system
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.datamodel import Particles

from amuse.community.sse.interface import SSE
from amuse.community.mesa.interface import MESA
from amuse.community.evtwin.interface import EVtwin
from amuse.community.evtwin2sse.interface import EVtwin2SSE

from amuse.community.huayno.interface import Huayno
from amuse.community.ph4.interface import ph4

from amuse.io import write_set_to_file

from amuse.couple.parallel_stellar_evolution import ParallelStellarEvolution
from amuse.couple.fallback_stellar_evolution import FallbackStellarEvolution

class notsorandom(object):
  def random(self,N):
    return numpy.array(range(N))/(N-1.)

def new_cluster(number_of_stars = 1000,Rscale=1. | units.parsec, seed=1234567):
    numpy.random.seed(seed)
    masses = new_salpeter_mass_distribution(
        number_of_stars, 
        mass_min = 0.3 | units.MSun,
        mass_max = 25.0 | units.MSun, 
        alpha = -2.35,random=notsorandom()
    )
    nbody_converter = nbody_system.nbody_to_si(masses.sum(), Rscale)
    particles = new_plummer_model(number_of_stars, nbody_converter)
    particles.mass = masses
    particles.radius=0. | units.RSun
    particles.move_to_center()
    cm,rcore,rhocore=particles.densitycentre_coreradius_coredens( unit_converter=nbody_converter )
    core=particles.cluster_core(unit_converter=nbody_converter,density_weighting_power=1)
    print core.radius.in_(units.parsec)
    print rcore.in_(units.parsec)
    print masses.max()
    print len(particles.bound_subset(unit_converter=nbody_converter))
    print particles.mass_segregation_Gini_coefficient(unit_converter=nbody_converter)
    return nbody_converter,particles
    
def evolve_with_stellar_evolution( grav,se,tend,conv,data,timestep=0.25| units.Myr,E0=None,dE=None):
    tnow=grav.model_time
    from_se_to_grav=se.particles.new_channel_to(grav.particles)
    if E0 is None:
      E0=grav.particles.kinetic_energy()+grav.particles.potential_energy()
    if dE is None:
      dE=0.*E0  
    while tnow<tend-timestep/2:
      grav.evolve_model(tnow+timestep)
      tnow=grav.model_time
      se.evolve_model(tnow)

      print tnow.in_(units.Myr),
      oldmass=grav.particles.mass
      from_se_to_grav.copy_attributes(["mass"])
      dmass=oldmass-grav.particles.mass
      parts=grav.particles.copy()
      de=dmass*(parts.specific_kinetic_energy()+parts.potential())
      
      cm_vel=parts.center_of_mass_velocity()
      grav.particles.velocity=grav.particles.velocity-cm_vel
      dE+=de.sum()+parts.total_mass()*(cm_vel[0]**2+cm_vel[1]**2+cm_vel[2]**2)/2

      E=parts.kinetic_energy()+parts.potential_energy()
      core=parts.cluster_core( unit_converter=conv,density_weighting_power=1 )
      
      print (E0-(E+dE))/E0,core.radius.in_(units.parsec),dmass.sum().in_(units.MSun),
      
      bs=parts.bound_subset(unit_converter=conv,density_weighting_power=1,core=core)
      gini=bs.mass_segregation_Gini_coefficient(unit_converter=conv,density_weighting_power=1,core=core)
      
      data.setdefault("bound mass",[]).append(bs.total_mass())
      data.setdefault("coreradius",[]).append(core.radius)
      data.setdefault("gini",[]).append(gini)
      data.setdefault("time",[]).append(tnow)
      
      print len(bs),
      print gini


      
class dummy_se(object):
  def __init__(self):
    self.particles=Particles()
  def evolve_model(self,tend):
    pass        

if __name__=="__main__":
  N=1000
  tend= 100. | units.Myr
  label="sse"


  def grav_code(conv):
    code=Huayno(conv, mode="openmp")
    code.parameters.epsilon_squared=(0.0 | units.parsec)**2
    code.timestep_parameter=0.001
    return code
  
  def fse():
    code=FallbackStellarEvolution(enforce_monotonic_mass_evolution=True, rms_weights=[5.,1.,1.])
    code._main_se.parameters.min_timestep_stop_condition=1. | units.yr
    code._main_se.parameters.maximum_number_of_stars=N/4
    return code

  def se_code():
#    code=ParallelStellarEvolution(fse,4)
    code=SSE()
    return code

  conv,particles=new_cluster(N,Rscale=4. | units.parsec,seed=123457)
  
  print conv.to_nbody(tend)

  grav=grav_code(conv)
  grav.particles.add_particles(particles)
  se=se_code()  
  se.particles.add_particles(particles)
  se.commit_particles()

  dump_times=1./16.*numpy.array(range(1601)) | units.Myr
  
  data=dict()
  core=grav.particles.cluster_core( unit_converter=conv )
  bs=grav.particles.bound_subset(unit_converter=conv,density_weighting_power=1,core=core)
  gini=bs.mass_segregation_Gini_coefficient(unit_converter=conv,density_weighting_power=1,core=core)
      
  data.setdefault("bound mass",[]).append(bs.total_mass())
  data.setdefault("coreradius",[]).append(core.radius)
  data.setdefault("gini",[]).append(gini)
  data.setdefault("time",[]).append(0.*tend)
  
  for tend in dump_times:
    evolve_with_stellar_evolution(grav,se,tend,conv,data,timestep=0.0625| units.Myr)

    i=int(tend/ (0.0625|units.Myr))
    write_set_to_file(grav.particles,label+"-grav-%6.6i"%i,"amuse",append_to_file=False)
    write_set_to_file(se.particles.copy(),label+"-se-%6.6i"%i,"amuse",append_to_file=False)
  
    f=open(label+"-data","wb")
    cPickle.dump(data,f)
    f.close()
