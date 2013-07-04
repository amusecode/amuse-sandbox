"""
  example of molecular cloud evolution with explictly 
  split SPH and grav evolution

  Initial condition is a smooth spherical cloud with random velocities
  as in Bonnell et al. (2003)  
  
"""  

import os
import numpy
import cPickle  
  
from matplotlib import pyplot 

from amuse.units import nbody_system
from amuse.units import units,constants
from amuse.datamodel import Particles

from amuse.community.fi.interface import Fi
from amuse.community.octgrav.interface import Octgrav
from fisink import SinkFi

from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from amuse.ext.derived_grav_systems import copycat
from amuse.ext.bridge import bridge

from amuse.io import write_set_to_file

from thermal_model import SPH_with_Thermal_Model, attach_thermal_model, SimplifiedThermalModel

import time

numpy.random.seed(123457)

snapshot_directory="./snapshots"

try:
  os.mkdir(snapshot_directory)
except:
  print snapshot_directory,"already exists"  


def generate_cloud(N=5000, Mcloud=10000. | units.MSun, Rcloud=1. | units.parsec,
           rhocloud=None,ndenscrit=64, ethep_ratio=0.01,power=-4,nf=64,ekep_ratio=1.):

    result=dict()

    if rhocloud is None:
      rhocloud=Mcloud/(4./3.*numpy.pi*Rcloud**3)
    else:
      Rcloud=((Mcloud/(4./3.*numpy.pi*rhocloud))**(1./3.)).in_(units.parsec)
    print "MCloud:", Mcloud
    print "RCloud:", Rcloud

    tff_cloud=(1/(4*numpy.pi*constants.G*rhocloud)**0.5).in_(units.Myr)

    result['t_ff']=tff_cloud

    conv = nbody_system.nbody_to_si(Mcloud,Rcloud)
 
    result['converter']=conv
 
    parts=molecular_cloud(targetN=N,nf=nf,convert_nbody=conv,ekep_ratio=ekep_ratio,
            base_grid=body_centered_grid_unit_cube,power=power,ethep_ratio=ethep_ratio).result

    result['particles']=parts

    rhomax=((1./6*numpy.pi**(5./2)*parts.u[0]**(3./2.)/constants.G**(3./2))/(parts.mass[0]*ndenscrit))**2

    rhomax=rhomax.in_(units.MSun/units.parsec**3)
    rhocloud=rhocloud.in_(units.MSun/units.parsec**3)
 
    result['rhomax']=rhomax 
    result['rhocloud']=rhocloud


    tff_rhomax=(1/(4*numpy.pi*constants.G*rhomax)**0.5).in_(units.Myr)
    torbit=((Rcloud**3/(constants.G*Mcloud))**0.5).in_(units.Myr)

    print "rho cloud:",rhocloud.in_(units.MSun/units.parsec**3),rhocloud.in_(units.g/units.cm**3),rhocloud.in_(units.amu/units.cm**3)
    print "resolution element:", parts.mass[0].in_(units.MSun),(ndenscrit*parts.mass[0]).in_(units.MSun)
    print "rhomax:",rhomax,rhomax.in_(units.amu/units.cm**3),rhomax.in_(units.g/units.cm**3)
    print "rhomax/rhocloud:",rhomax/rhocloud
    print
    print "t_ff (rhomax):",tff_rhomax
    print "t_ff (rhocloud):",tff_cloud
    print "t_orbit:", torbit
    print
    
    meanu=parts.u.mean()
    csound=(meanu**0.5).in_(units.km/units.s)
    temperature=((2.2| units.amu)*meanu/constants.kB).in_(units.K)
    print "soundspeed:", csound
    print "temperature:", temperature

    result['temperature']=temperature

    return result


    
def run_mc(conv, thermal_model, istart, snapshot_times_tff, dt, rhomax, parts=None,
                      merge_radius=100. | units.AU,bridge_selfgrav=False, isothermal=True,dryrun=False):

    print "dt:", dt.in_(units.Myr)
    tend=snapshot_times[-1].in_(units.Myr)
    print "tend:",tend


    sph=attach_thermal_model(SinkFi)(conv,mode="openmp",redirection="none",
               density_threshold=rhomax,merge_radius=merge_radius, 
               thermal_model=thermal_model, timestep=dt/2, calc_net_luminosity=True)

    sph.parameters.use_hydro_flag=True
    sph.parameters.radiation_flag=False
    sph.parameters.self_gravity_flag=False if bridge_selfgrav else True
    if isothermal:
      sph.parameters.gamma=1
      sph.parameters.isothermal_flag=True
      sph.parameters.integrate_entropy_flag=False
    else:
      sph.parameters.gamma=1.4
      sph.parameters.isothermal_flag=False
      sph.parameters.integrate_entropy_flag=True      
    sph.parameters.timestep=dt/2
    sph.parameters.verbosity = 0

    sinks=Particles()
    if parts is None or istart>=0:
      parts=read_set_from_file(snapshot_directory+'/gas-%6.6i'%istart,'amuse')
      try:
        sinks=read_set_from_file(snapshot_directory+'/sink-%6.6i'%istart,'amuse')
      except:
        pass

    sph.gas_particles.add_particles(parts)
    if len(sinks)>0:
      sph.dm_particles.add_particles(sinks)
      sph.reinit_sink(sinks)

    print
    print "smoothing:", sph.gas_particles.h_smooth.mean().in_(units.AU)
    rhoh3=(sph.gas_particles.density*sph.gas_particles.h_smooth**3).mean()
    print "minimum smoothing:", ((rhoh3/rhomax)**(1./3)).in_(units.AU)
    print
    
    if bridge_selfgrav:
      print "bridge selfgrav"
      def Fi_mp(conv):
        return Fi(conv,mode="openmp",redirection="none")
    
      grav=copycat(Fi_mp, sph, conv)

      sys=bridge(verbose=False)
      sys.add_system(sph,(grav,),False)
      sys.timestep=dt
    else:
      sys=sph

    if istart<0: 
      data=dict()
      data['time']=[]
      data['luminosity']=[]
      data['radiated_energy']=[]
    else:
      f=open(snapshot_directory+"/data.pkl","rb")
      cPickle.load(f)
      f.close()

    i=istart
    t0=time.time()

    if dryrun:
      return

    for ttarget in snapshot_times[istart+1:]:
        i=i+1
        sys.evolve_model(ttarget)

        gas=sph.gas_particles.copy()
        sinks=sph.sink_particles.copy()
        
        gas.collection_attributes.time=sys.model_time
        sinks.collection_attributes.time=sys.model_time
        
        write_set_to_file(gas,snapshot_directory+'/gas-%6.6i'%i,'amuse',append_to_file=False)
        write_set_to_file(sinks,snapshot_directory+'/sink-%6.6i'%i,'amuse',append_to_file=False)
        
        data['time'].append(sys.model_time)
        data['luminosity'].append(sph.total_luminosity)
        data['radiated_energy'].append(sph.radiated_energy)
        
        f=open(snapshot_directory+"/data.pkl","wb")
        cPickle.dump(data,f)
        f.close()
        
        print ttarget.in_(units.Myr),
        print sph.gas_particles.rho.amax().in_(units.MSun/units.parsec**3),
        print len(sph.dm_particles)
        
        t2=time.time()
        
        if i>0:
          eta=(t2-t0)/ttarget*(tend-ttarget)
          print "ETA (hours):", eta/3600

  
if __name__ in ("__main__","__plot__"):


    istart=-1
#    nsnap=12
#    snapshot_times_tff=3.*numpy.arange(nsnap+1)/nsnap
    snapshot_times_tff=[0., 0.5, 1.0, 1.5, 2.0, 2.5]

   
#generate_cloud(N=5000, Mcloud=10000. | units.MSun, Rcloud=1. | units.parsec,
#           rhocloud=None, rhomax=None,ndenscrit=64, 
#           dt_tff=0.01,merge_radius=100. | units.AU,
#           ethep_ratio=0.01,power=-4,nf=64,ekep_ratio=1.)
   
    ic=generate_cloud(N=100000, Mcloud=1000. | units.MSun, Rcloud=1. | units.parsec,
             ethep_ratio=0.02)
        
    snapshot_times=snapshot_times_tff*ic['t_ff']
    dt=0.01*ic['t_ff']
    parts=ic['particles']
    conv=ic['converter']
    rhomax=ic['rhomax']
    Tmin=ic['temperature']
    
    def thermal_model():
      return SimplifiedThermalModel(Tmin=Tmin)
        
#run_mc(conv, thermal_model, istart, snapshot_times_tff, dt, rhomax, parts=None,
#                      merge_radius=100. | units.AU,bridge_selfgrav=False, isothermal=True,dryrun=False):

    run_mc(conv, thermal_model, istart, snapshot_times, dt, rhomax, parts=parts,dryrun=False,bridge_selfgrav=True)

