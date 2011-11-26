import os
import sys
import numpy
import time

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from grav_gas_sse import grav_gas_sse

from amuse.units import nbody_system
from amuse.units import units
from amuse.units import constants

from amuse.community.fi.interface import Fi
from boxedfi import BoxedFi
from amuse.community.octgrav.interface import Octgrav
from amuse.community.gadget2.interface import Gadget2
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPE

from SSEplus import SSEplus

from amuse import datamodel

from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.plummer import new_plummer_sphere
from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.evrard_test import body_centered_grid_unit_cube

import cProfile

#numpy.random.seed(12345)

def smaller_nbody_power_of_two(dt, conv):
  nbdt=conv.to_nbody(dt).value_in(nbody_system.time)
  idt=numpy.floor(numpy.log2(nbdt))
  return conv.to_si( 2**idt | nbody_system.time)


def clustergas_restart(runid, snapshot, t_end=30. | units.Myr,
                         dt_plot= 0.05 | units.Myr,newid=None,
                         new_gas_options=(),
                 grav_code_extra=dict(mode='gpu', redirection='none'),
                 gas_code_extra=dict(number_of_workers=3,use_gl=False, redirection='none'),
                 se_code_extra=dict(redirection='none'),
                 grav_couple_code_extra=dict()):
  if runid is None:
    raise Exception                       
  print "restarting run "+runid+" at snapshot %i"%snapshot                       
  print
  
  conv,sys=grav_gas_sse.load_system_state(runid+"/dump-%6.6i" %snapshot,
             new_gas_options,
             grav_code_extra=grav_code_extra,
             gas_code_extra=gas_code_extra,
             se_code_extra=se_code_extra,
             grav_couple_code_extra=grav_couple_code_extra)

  dt =smaller_nbody_power_of_two(dt_plot, conv)

  if newid is not None:
    runid=newid
    print "switch to:",runid
    try:
      os.mkdir(runid)
    except:
      raise Exception

  print
  print "so far so good...resuming evolve!"

  nsnap=snapshot
  t=sys.model_time
  
  while (t < t_end-dt/2):
    sys.evolve_model(t+dt)
    sys.synchronize_model()
    t=sys.model_time
    tout=(t).value_in(units.Myr)
    ek=sys.kinetic_energy.value_in(1.e51*units.erg)
    ep=sys.potential_energy.value_in(1.e51*units.erg)
    eth=sys.thermal_energy.value_in(1.e51*units.erg)
    ef=sys.feedback_energy.value_in(1.e51*units.erg)
    print 't Ek Ep Eth Ef:', tout,ek,ep,eth,ef,ek+ep+eth-ef
    nsnap+=1
    sys.dump_system_state(runid+"/dump-%6.6i" %nsnap)
                         

#original feedback_efficiency=0.01
def clustergas(sfeff=0.05,Nstar=1000,Ngas=1000, t_end=30. | units.Myr,
                 dt_plot= 0.05 | units.Myr, Rscale= 0.5 | units.parsec,
                 runid="runtest",feedback_efficiency=0.03, subvirialfac=1.,
                 grav_code=PhiGRAPE,gas_code=Gadget2,se_code=SSEplus,grav_couple_code=Octgrav,
                 grav_code_extra=dict(mode='gpu', redirection='none'),
                 gas_code_extra=dict(number_of_workers=3,use_gl=False, redirection='none'),
                 se_code_extra=dict(redirection='none'),
                 grav_couple_code_extra=dict()
                 ):

  try:
    os.mkdir(runid)
  except:
    pass  
  print "Nstar, Ngas:", Nstar,Ngas
  eps=0.05 * Rscale
  eps_star=0.001 * Rscale

  star_masses = new_salpeter_mass_distribution(Nstar, mass_max = 100. | units.MSun)
  total_star_mass = star_masses.sum()
  total_mass=total_star_mass/sfeff
  print "maxmass:", max(star_masses)
  print "Tcross", (2.8*(Rscale**3/ constants.G/ total_star_mass)**0.5).in_(units.Myr)

  print sorted(star_masses)[-10:]
#  raise Exception

  conv = nbody_system.nbody_to_si(total_mass,Rscale)

  star_parts=new_plummer_sphere(Nstar,convert_nbody=conv)
  star_parts.mass=star_masses
  star_parts.radius=eps_star
  
# sub-virialized
  star_parts.vx=star_parts.vx*subvirialfac
  star_parts.vy=star_parts.vy*subvirialfac
  star_parts.vz=star_parts.vz*subvirialfac

  print "total cluster mass:", total_mass.in_(units.MSun)
  print "star mass:", total_star_mass.in_(units.MSun)
  print "gas mass:", (total_mass-total_star_mass).in_(units.MSun)

  print "t_end:", conv.to_nbody(t_end)
  
  gas_parts=new_plummer_gas_model(Ngas,convert_nbody=conv, base_grid=body_centered_grid_unit_cube)
  gas_parts.h_smooth=0. | units.parsec
  gas_parts.mass=gas_parts.mass*(1-sfeff)

  print "gas particle mass:",  ((total_mass-total_star_mass)/len(gas_parts)).in_(units.MSun)

  mu=1.4 | units.amu
  gamma1=1.6667-1
#  print 'min Temp:', (gamma1*min(gas_parts.u)*(1.4*units.amu)/constants.kB).in_(units.K)
  print 'min Temp:', (gamma1*min(gas_parts.u)*mu/constants.kB).in_(units.K)

  mgas=(total_mass-total_star_mass)/len(gas_parts)
  print max(gas_parts.u)**0.5
  
  dt =smaller_nbody_power_of_two(dt_plot, conv)
  dt_star=smaller_nbody_power_of_two( 0.001 | units.Myr, conv) #0.001

#original settings  
#  dt_sph=dt_star
#  dt_fast=dt_star*8  # 8
#  dt_feedback=dt_fast*2 # 1

#new, improved, more snapshots (thanks inti)
  dt_fast=dt / 8  # 8
  dt_feedback=dt # 1

  feedback_dt=(10. | units.yr, dt_fast)

  if not dt_star<=dt_fast<=dt_feedback:
    raise Exception 
     
  print 'dt_plot:', conv.to_nbody(dt_plot), dt_plot.in_(units.Myr)
  print 'dt:', conv.to_nbody(dt), dt.in_(units.Myr)
  print 'dt_feedback:', conv.to_nbody(dt_feedback), dt_feedback.in_(units.Myr)
  print 'dt_fast:', conv.to_nbody(dt_fast), dt_fast.in_(units.Myr)
  print 'dt_star:', conv.to_nbody(dt_star), dt_star.in_(units.Myr)

  star_parts.move_to_center()
  gas_parts.move_to_center()
  
  sys=grav_gas_sse(grav_code,gas_code,se_code,grav_couple_code,
               conv,mgas,star_parts,gas_parts,dt_feedback,dt_fast,
               grav_parameters=(["epsilon_squared", eps_star**2],["timestep_parameter", 0.001]),
#                                ["timestep", dt_star]),
               gas_parameters=(["time_max", 32. | units.Myr],
#                               ["courant", 0.15 | units.none],
                               ["n_smooth", 64 | units.none],
#                               ["artificial_viscosity_alpha",1.| units.none],
                               ["n_smooth_tol", 0.005 |units.none],
                               ## NB
#                               ["max_size_timestep",7500. | units.yr],
                               ## NB
                               ["time_limit_cpu", 3600000 | units.s]),
               couple_parameters=(["epsilon_squared", eps**2],
                            ["opening_angle",0.5]),
               feedback_efficiency=feedback_efficiency, feedback_radius=0.025*Rscale,
               feedback_dt=feedback_dt,
               grav_code_extra=grav_code_extra,
               gas_code_extra=gas_code_extra,
               se_code_extra=se_code_extra,
               grav_couple_code_extra=grav_couple_code_extra
               )

  nsnap=0
  sys.synchronize_model()
  t=sys.model_time
  tout=t.value_in(units.Myr)
  ek=sys.kinetic_energy.value_in(1.e51*units.erg)
  ep=sys.potential_energy.value_in(1.e51*units.erg)
  eth=sys.thermal_energy.value_in(1.e51*units.erg)
  ef=sys.feedback_energy.value_in(1.e51*units.erg)
  print 't Ek Ep Eth Ef:', tout,ek,ep,eth,ef,ek+ep+eth-ef
  sys.dump_system_state(runid+"/dump-%6.6i" %nsnap)

  print 'max smooth:', max(sys.sph.gas_particles.radius).in_(units.parsec)
  print 'mean smooth:', sys.sph.gas_particles.radius.mean().in_(units.parsec)
  print 'min smooth:', min(sys.sph.gas_particles.radius).in_(units.parsec)
  print 'eps:',eps.in_(units.parsec)
  print 'feedback radius:',(0.01*Rscale).in_(units.parsec)

  while (t<t_end-dt/2):
    beginning = time.time()
    sys.evolve_model(t+dt)
    sys.synchronize_model()
    t=sys.model_time
    tout=t.value_in(units.Myr)
    ek=sys.kinetic_energy.value_in(1.e51*units.erg)
    ep=sys.potential_energy.value_in(1.e51*units.erg)
    eth=sys.thermal_energy.value_in(1.e51*units.erg)
    ef=sys.feedback_energy.value_in(1.e51*units.erg)
    print 't Ek Ep Eth Ef:', tout,ek,ep,eth,ef,ek+ep+eth-ef
    nsnap+=1
    sys.dump_system_state(runid+"/dump-%6.6i" %nsnap)

    end = time.time()

    print 'iteration', nsnap , 'took:', (end - beginning), 'seconds'


def generate_cluster_with_gas(IMF,sfeff,Nstar,Rscale,mgas,
                              subvirialfac=1.,t_end=10| units.Myr):
  
    total_star_mass,star_mass = IMF.next_set(Nstar)

    total_mass=total_star_mass/sfeff
    total_gas_mass=(1-sfeff)*total_mass
    Ngas=long((total_gas_mass/mgas).value_in(units.none))

    print "Ngas:",Ngas  
    print "Tcross", (2.8*(Rscale**3/ constants.G/ total_mass)**0.5).in_(units.Myr)
    print sorted(star_mass)[-10:]
  
    conv = nbody_system.nbody_to_si(total_mass,Rscale)

    star_parts=MakePlummerModel(Nstar,convert_nbody=conv).result
    star_parts.mass=star_mass
  
# sub-virialized
    star_parts.vx=star_parts.vx*subvirialfac
    star_parts.vy=star_parts.vy*subvirialfac
    star_parts.vz=star_parts.vz*subvirialfac

    print "total cluster mass:", total_mass.in_(units.MSun)
    print "star mass:", total_star_mass.in_(units.MSun)
    print "gas mass:", (total_mass-total_star_mass).in_(units.MSun)

    print "t_end in nbody:", conv.to_nbody(t_end)
  
    
    gas_parts=MakePlummerGasModel(Ngas,convert_nbody=conv, base_grid=body_centered_grid_unit_cube).result
    print gas_parts.mass[0]*(1-sfeff),mgas
    gas_parts.mass=mgas

    star_parts.move_to_center()
    gas_parts.move_to_center()

    return star_parts,gas_parts


def merge_clusters_with_gas(sfeff=[0.3,0.3], Nstar= [1000,1000],
                 totalNgas=40000, t_end=30. | units.Myr,
                 dt_plot= 0.05 | units.Myr, Rscale= units.parsec([0.25,0.25]),
                 tmerge=2.5 | units.Myr,eps_star=0.001 | units.parsec,
                 runid="runtest",feedback_efficiency=0.01, subvirialfac=1.):

  try:
    os.mkdir(runid)
  except:
    pass  
  print "Nstar, Ngas:", Nstar,totalNgas
  
  eps=0.05 * Rscale.mean()

  IMF=SalpeterIMF(mass_max = 100. | units.MSun)
  
  mean_star_mass=IMF.mass_mean()

  approx_gas_mass=0. * mean_star_mass 
  for e,N in zip(sfeff,Nstar):
    approx_gas_mass+=(1-e)*N*mean_star_mass/e
  mgas=approx_gas_mass/totalNgas

  print "gas particle mass:", mgas.in_(units.MSun)
  
  stars=[]
  gas=[]
  total_mass=0. | units.MSun
  total_gas_mass=0. | units.MSun
  total_star_mass=0. | units.MSun
  for e,Ns,r in zip(sfeff,Nstar,Rscale):
   s,g=generate_cluster_with_gas(IMF,e,Ns,r,mgas,
                                         subvirialfac=subvirialfac,t_end=t_end)
   stars.append(s)
   gas.append(g)
   total_mass+=s.mass.sum()+g.mass.sum()
   total_star_mass+=s.mass.sum()
   total_gas_mass+=g.mass.sum()

  print "total mass:", total_mass.in_(units.MSun)
  print "total star mass:", total_star_mass.in_(units.MSun)
  print "total gas mass:", total_gas_mass.in_(units.MSun)

  print "t_end:", t_end

  Rsep=(constants.G*total_mass*tmerge**2)**(1./3)

  stars[0].x=stars[0].x-Rsep/2
  stars[1].x=stars[1].x+Rsep/2
  gas[0].x=gas[0].x-Rsep/2
  gas[1].x=gas[1].x+Rsep/2

  conv = nbody_system.nbody_to_si(total_mass,Rsep)

  print "tmerge:",tmerge
  print "Rsep:",Rsep.in_(units.parsec)

  allgas=datamodel.Particles(0)
  allstars=datamodel.Particles(0)
  allgas.add_particles(gas[0])
  allgas.add_particles(gas[1])
  allstars.add_particles(stars[0])
  allstars.add_particles(stars[1])
  allgas.h_smooth=0.| units.parsec
  allstars.radius=eps_star


  mu=1.4 | units.amu
  gamma1=1.6667-1
#  print 'min Temp:', (gamma1*min(allgas.u)*(1.4*units.amu)/constants.kB).in_(units.K)
  print 'min Temp:', (gamma1*min(allgas.u)*mu/constants.kB).in_(units.K)

  print max(allgas.u)**0.5
  
  dt =smaller_nbody_power_of_two(dt_plot, conv)
  dt_star=smaller_nbody_power_of_two( 0.001 | units.Myr, conv) #0.001
  dt_sph=dt_star
  dt_fast=dt_star  # 8
  dt_feedback=dt_fast*4 # 2
    
  if not dt_star<=dt_fast<=dt_feedback:
    raise Exception 
     
  print 'dt_plot:', conv.to_nbody(dt_plot), dt_plot.in_(units.Myr)
  print 'dt:', conv.to_nbody(dt), dt.in_(units.Myr)
  print 'dt_feedback:', conv.to_nbody(dt_feedback), dt_feedback.in_(units.Myr)
  print 'dt_fast:', conv.to_nbody(dt_fast), dt_fast.in_(units.Myr)
  print 'dt_star:', conv.to_nbody(dt_star), dt_star.in_(units.Myr)
  print 'dt_sph:', conv.to_nbody(dt_sph), dt_sph.in_(units.Myr)
  
  sys=grav_gas_sse(PhiGRAPE,Gadget2,SSEplus,Octgrav, 
               conv,mgas,allstars,allgas,dt_feedback,dt_fast,
               grav_parameters=(["epsilon_squared", eps_star**2],["timestep_parameter", 0.001]),
#                                ["timestep", dt_star]),
               gas_parameters=(["time_max", 32. | units.Myr],
#                               ["courant", 0.15 | units.none],
                               ["n_smooth", 64 | units.none],
#                               ["artificial_viscosity_alpha",1.| units.none],
                               ["n_smooth_tol", 0.005 |units.none],
                               ## NB
#                               ["max_size_timestep",7500. | units.yr],
                               ## NB
                               ["time_limit_cpu", 3600000 | units.s]),
               couple_parameters=(["epsilon_squared", eps**2],
                            ["opening_angle",0.5]),
               feedback_efficiency=feedback_efficiency, feedback_radius=0.025*Rscale)

  nsnap=0
  sys.synchronize_model()
  t=sys.model_time
  tout=t.value_in(units.Myr)
  ek=sys.kinetic_energy.value_in(1.e51*units.erg)
  ep=sys.potential_energy.value_in(1.e51*units.erg)
  eth=sys.thermal_energy.value_in(1.e51*units.erg)
  ef=sys.feedback_energy.value_in(1.e51*units.erg)
  print 't Ek Ep Eth Ef:', tout,ek,ep,eth,ef,ek+ep+eth-ef
  sys.dump_system_state(runid+"/dump-%6.6i" %nsnap)

  print 'max smooth:', max(sys.sph.gas_particles.radius).in_(units.parsec)
  print 'mean smooth:', sys.sph.gas_particles.radius.mean().in_(units.parsec)
  print 'min smooth:', min(sys.sph.gas_particles.radius).in_(units.parsec)
  print 'eps:',eps.in_(units.parsec)
  print 'feedback radius:',(0.01*Rscale).in_(units.parsec)

  while (t<t_end-dt/2):
    sys.evolve_model(t+dt)
    sys.synchronize_model()
    t=sys.model_time
    tout=t.value_in(units.Myr)
    ek=sys.kinetic_energy.value_in(1.e51*units.erg)
    ep=sys.potential_energy.value_in(1.e51*units.erg)
    eth=sys.thermal_energy.value_in(1.e51*units.erg)
    ef=sys.feedback_energy.value_in(1.e51*units.erg)
    print 't Ek Ep Eth Ef:', tout,ek,ep,eth,ef,ek+ep+eth-ef
    nsnap+=1
    sys.dump_system_state(runid+"/dump-%6.6i" %nsnap)

  
if __name__=="__main__":
  import time
#  cProfile.run("clustergas()","prof")
  t1=time.time()
  nb=merge_clusters_with_gas()
  t2=time.time()
  print 'time:',  t2-t1


