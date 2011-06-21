import numpy
import cPickle

from matplotlib import pyplot

from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.community.fi.interface import Fi
from amuse.support.data.core import Grid


from amuse.support.io import read_set_from_file

first=350
nsnap=351
cm=None
N=100

def fig(i):
  global cm
  g=read_set_from_file('dump-%6.6i'%i,'amuse')

  t=i*1.
  
  conv=nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)

  fi=Fi(convert_nbody=conv,use_gl=False,mode='periodic')

  fi.parameters.use_hydro_flag=True
  fi.parameters.radiation_flag=False
  fi.parameters.self_gravity_flag=False
  fi.parameters.gamma=1
  fi.parameters.isothermal_flag=True
  fi.parameters.integrate_entropy_flag=False
  fi.parameters.timestep=.01 | units.Myr  
  fi.parameters.verbosity=0
  fi.parameters.periodic_box_size=30 | units.kpc
  fi.commit_parameters()

  fi.gas_particles.add_particles(g)
  fi.commit_particles()

  L=30.
  x,y=numpy.indices( ( N+1,N+1 ))

  x=L*(x.flatten()-N/2.)/N
  y=L*(y.flatten()-N/2.)/N
  z=x*0.
  vx=0.*x
  vy=0.*x
  vz=0.*x

  x=units.kpc(x)
  y=units.kpc(y)
  z=units.kpc(z)
  vx= (units.km/units.s).new_quantity(vx)
  vy= (units.km/units.s).new_quantity(vy)
  vz= (units.km/units.s).new_quantity(vz)

  rho,rhovx,rhovy,rhovz,rhoe=fi.get_hydro_state_at_point(y,x,z,vx,vy,vz)

  rho=rho.reshape((N+1,N+1))
  rhoe=rhoe.reshape((N+1,N+1))

#  print ((rho.value_in(units.amu/units.cm**3))).min()
#  print ((rho.value_in(units.amu/units.cm**3))).max()
#  print (numpy.log(rho.value_in(1.67e-24*units.g/units.cm**3))).min()
#  print (numpy.log(rho.value_in(1.67e-24*units.g/units.cm**3))).max()
  f=pyplot.figure(figsize=(6,6))
  pyplot.imshow(numpy.log10(1.e-20+rho.value_in(units.amu/units.cm**3)),
    extent=[-30.,30.,-30.,30],origin='lower',vmin=-5,vmax=-2,cmap='copper')
  pyplot.figtext(0.15,0.04,t)
  pyplot.savefig('dens-%6.6i.eps'%i)
  f.clear()
  pyplot.close(f)


if __name__=='__main__':
  print
  for i in range(first,nsnap):
    print '.',
    fig(i)
