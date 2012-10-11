import os
import sys
import numpy

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    
from amuse.units import nbody_system
from amuse.units import units
from amuse.io import write_set_to_file    
    
from amuse.community.fi.interface import Fi
from amuse.community.huayno.interface import Huayno
from amuse.community.ph4.interface import ph4
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE


from amuse.ic.fractalcluster import new_fractal_cluster_model

from amuse.ext.gasplummer import MakePlummerGasModel
from amuse.ic.plummer import MakePlummerModel
from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.evrard_test import body_centered_grid_unit_cube
from amuse.ext.radial_profile import radial_density

from fast import FAST
from copycat import reinitializecopycat
 
import logging
#logging.basicConfig(level=logging.DEBUG)

def energy_plot(time,ek,ep,eth,label="none"):
  if not HAS_MATPLOTLIB:
    return
    
  f=pyplot.figure(figsize = (8, 6))
  pyplot.xlabel(r'time')
  pyplot.ylabel(r'energy')
  pyplot.plot(time,ek,'g:')
  pyplot.plot(time,ep,'y--')
  pyplot.plot(time,eth,'b-.')
  pyplot.plot(time,ek+ep+eth,'r')
  pyplot.savefig(label+"-ekepeth.eps")
  pyplot.close(f)


def radial_plot(time,rad_dens,label="none"):
  c=['g','b','y','c','r:','g:']

  f=pyplot.figure(figsize = (8, 6))
  pyplot.xlabel(r'radius')
  pyplot.ylabel(r'density')
  for i in range(len(rad_dens)):
    pyplot.loglog(rad_dens[i][0],rad_dens[i][1],c[i])
  
  ascl=1/1.695
  ra,dens=rad_dens[0]
  pyplot.loglog(ra, 3./4./numpy.pi/ascl**3/(1+(ra**2/ascl**2))**(5./2),'r')
  pyplot.xlim(0.1,10)
  pyplot.ylim(1.e-7,10)
  pyplot.savefig(label+'-rad_dens.eps')
  pyplot.close(f)
  
  
def mixed_fractal_model(N,Ngas,fgas=0.7,fd=1.6,ethep_ratio=0.01):
  
  gas = new_fractal_cluster_model(N=Ngas,fractal_dimension=fd,random_seed=12345321,virial_ratio=0.25)
  stars = new_fractal_cluster_model(N=N,fractal_dimension=fd,random_seed=12345321,virial_ratio=0.25)
  gas.hsmooth=0.01 | nbody_system.length
  stars.radius=0.01 | nbody_system.length

#  ep=gas.potential_energy(G=nbody_system.G)
#  gas.u=ethep_ratio*ep/gas.total_mass()
  gas.u=0.001 | nbody_system.specific_energy
  for g in gas:
    g.position+=0.05*nbody_system.length(numpy.random.random(3)-0.5)   
  stars.mass*=(1-fgas)
  gas.mass*=fgas
  
  return stars,gas
  
def xyplot(gas,stars,i):
  f=pyplot.figure(figsize=(8,8))
  pyplot.plot(gas.x.number,gas.y.number,'g.')
  pyplot.plot(stars.x.number,stars.y.number,'r+')
  pyplot.xlim(-1.,1.)
  pyplot.ylim(-1.,1.)
  pyplot.savefig("stargas_xy-%6.6i.png"%i)
  f.clear()
  pyplot.close()
  
def run_mixed_fractal(N=100,Ngas=10000,fgas=0.7,fd=1.6,label="none"):

  dt=1./32 | nbody_system.time
  dt_fast=1./128 | nbody_system.time
  tend=4. | nbody_system.time

  stars,gas=mixed_fractal_model(N,Ngas,fgas=fgas,fd=fd)
    
  sph=Fi()

  sph.parameters.use_hydro_flag=True
  sph.parameters.radiation_flag=False
  sph.parameters.self_gravity_flag=False
  sph.parameters.integrate_entropy_flag=True
  sph.parameters.timestep=dt_fast/2  
  sph.parameters.verbosity=99

  sph.commit_parameters()

  sph.gas_particles.add_particles(gas)

  grav=Hermite()
  
  grav.parameters.epsilon_squared=(0.01 | nbody_system.length)**2
  
  grav.particles.add_particles(stars)
   
  couple1=reinitializecopycat(Fi, (sph,),None,parameters=(["epsilon_squared",(0.01 | nbody_system.length)**2],),extra={"channel_type": "sockets"})
  
  fast=FAST(verbose=True)
  fast.timestep=dt_fast
  fast.add_system(sph,(grav,couple1))
  fast.add_system(grav,(couple1,))
    
  tstart=fast.model_time
   
  fast.synchronize_model()
  time,Ek,Ep,Eth,rad_dens=[],[],[],[],[]
  time.append(tstart.number)
  ek=fast.kinetic_energy
  Ek.append(ek.number)
  ep=fast.potential_energy
  Ep.append(ep.number)
  eth=fast.thermal_energy
  Eth.append(eth.number)
  
  r=(fast.gas_particles.x**2+fast.gas_particles.y**2+fast.gas_particles.z**2)**0.5
  ra,dens=radial_density(r,fast.gas_particles.mass,500)
  ra=numpy.array(map(lambda x: x.number,ra))
  dens=numpy.array(map(lambda x: x.number,dens))
  rad_dens.append([ra,dens])

  xyplot(sph.gas_particles,grav.particles,0)

  tnow=tstart
  i=0
  while tnow<tend-dt/2:
    i+=1
    print tnow
    tnow=tnow+dt
    fast.evolve_model(tnow)
    fast.synchronize_model()
    tcur=fast.model_time
    time.append(tnow.number)
    ek=fast.kinetic_energy
    Ek.append(ek.number)
    ep=fast.potential_energy
    Ep.append(ep.number)
    eth=fast.thermal_energy
    Eth.append(eth.number)

    r=(fast.gas_particles.x**2+fast.gas_particles.y**2+fast.gas_particles.z**2)**0.5
    ra,dens=radial_density(r,fast.gas_particles.mass,500)
    ra=numpy.array(map(lambda x: x.number,ra))
    dens=numpy.array(map(lambda x: x.number,dens))
    rad_dens.append([ra,dens])
    
    p=sph.gas_particles.copy()
    write_set_to_file(p,'gas-%6.6i'%i)
    p=grav.particles.copy()
    write_set_to_file(p,'stars-%6.6i'%i)

    xyplot(sph.gas_particles,grav.particles,i)

  time=numpy.array(time)
  Ek=numpy.array(Ek)
  Ep=numpy.array(Ep)
  Eth=numpy.array(Eth)

  energy_plot(time,Ek,Ep,Eth,label=label)
  radial_plot(time,rad_dens,label=label)

if __name__=="__main__":
  run_mixed_fractal(N=100,Ngas=10000,fgas=0.7,label="mixed_fractal-10k")
