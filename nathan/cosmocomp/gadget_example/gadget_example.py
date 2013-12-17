import numpy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

from amuse.units import generic_unit_converter
from amuse.units import units

from amuse.community.gadget2.interface import Gadget2

from box import box

def make_map(sph,L,N=100):

    x,y=numpy.indices( ( N+1,N+1 ))

    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=x*0.
        
    vx=units.kms(numpy.zeros_like(x.number))
    vy=units.kms(numpy.zeros_like(x.number))
    vz=units.kms(numpy.zeros_like(x.number))

    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    rho=rho.reshape((N+1,N+1))

    return numpy.transpose(rho)

if __name__=="__main__":

  N=10000
  L= 10| units.parsec
  rho= 1.14*1 | units.amu/units.cm**3
  u=(5.e11 | units.erg/units.g).in_(units.cm**2/units.s**2) # =5000K
  
  tend=1. | units.Myr
  dt=10000 | units.yr
  
  print (u**0.5).in_(units.kms)
  print ((L/u**0.5)/dt)
    
  particles=box(N, L,rho,u)

  UnitLength=L
  UnitMass=particles.mass.sum()
  UnitVelocity=units.kms
  
  convert=generic_unit_converter.ConvertBetweenGenericAndSiUnits(UnitLength, UnitMass, UnitVelocity)
  sph=Gadget2(convert,mode='periodic_nogravity')#,redirection='none')
  
  sph.parameters.periodic_box_size=L
  
  sph.gas_particles.add_particles(particles)
  
  i=0
  t=0. | units.Myr
  while t<(tend-dt/2):
    t=t+dt
    i=i+1
    sph.evolve_model(t)
    print t.in_(units.Myr),sph.model_time.in_(units.Myr)
    rho=make_map(sph,L)
    f=pyplot.figure(figsize=(8,8))
    LL=L.number
    pyplot.imshow(numpy.log10(rho.value_in(units.amu/units.cm**3)),
        extent=[-LL/2,LL/2,-LL/2,LL/2],vmin=0,vmax=2,origin='lower')
    pyplot.savefig('map-%6.6i.png'%i)
    f.clear()
    pyplot.close(f)

    
