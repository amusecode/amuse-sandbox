import numpy

from math import *

from amuse.units import nbody_system
from amuse.units import units
from amuse.units.quantities import zero

from amuse.datamodel import Particles

from amuse.ext.evrard_test import *

def box(N, L, rho,u,base_grid=uniform_random_unit_cube):

  x,y,z=base_grid(N).make_xyz()
  Nresult=len(x)
  part=Particles(Nresult)
  part.x=L*x
  part.y=L*y
  part.z=L*z
  part.mass=(rho*L**3)/Nresult

  vunit=(u.unit)**0.5
  part.vx=vunit(numpy.zeros_like(x))
  part.vy=vunit(numpy.zeros_like(x))
  part.vz=vunit(numpy.zeros_like(x))
  
  part.u=u
  part.radius=L/Nresult**(1./3.)
  
  return part

if __name__=="__main__":

  N=100
  L=100| units.parsec
  rho= 1.14*1 | units.amu/units.cm**3
  u=(5.e11 | units.erg/units.g).in_(units.cm**2/units.s**2) # =5000K
  
  part=box(N, L,rho,u)
  print part
  print len(part)
