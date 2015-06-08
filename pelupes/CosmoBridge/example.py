from amuse.units import units,nbody_system

from amuse.ext.cosmo import Cosmology

from CosmoBridge import VacuumBoundary_CosmoBridge

from amuse.community.bhtree.interface import BHTree
from amuse.community.fi.interface import Fi

from amuse.io import read_set_from_file

from matplotlib import pyplot

if __name__=="__main__":
  
  z0=100.
  zcollaps=30.

  cosmo=Cosmology()

  print cosmo.agefromz(zcollaps).in_(units.Myr)
  print cosmo.agefromz(z0).in_(units.Myr)

  sys=read_set_from_file("example.amuse","amuse")
  
  conv=nbody_system.nbody_to_si(1.e6 |units.MSun,10. | units.parsec)
  
  grav=Fi(conv)
  
  grav.parameters.epsilon_squared=(5 | units.parsec)**2
  grav.particles.add_particles(sys)
  
  bridge=VacuumBoundary_CosmoBridge(z0, cosmo)
  bridge.add_system(grav)
  
  tend=cosmo.agefromz(zcollaps)
  dt=0.5 | units.Myr
  tnow=bridge.time
  
  pyplot.ion()
  pyplot.figure()
  pyplot.show()
  
  while bridge.time<tend:
      bridge.evolve_model(tnow+dt)
      tnow=bridge.time
      a=cosmo.afromage(tnow)
      z=1./a-1.
      print tnow.in_(units.Myr),z

      pyplot.clf()
      x=bridge.particles.x.value_in(units.parsec)
      y=bridge.particles.y.value_in(units.parsec)
      pyplot.plot(x,y,'k.')
      pyplot.draw()
      
