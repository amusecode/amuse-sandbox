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
from amuse.units import constants
from amuse.units.quantities import zero

from amuse.legacy.fi.interface import Fi

from amuse.ext.salpeter import SalpeterIMF
from amuse.ic.plummer import new_plummer_sphere
numpy.random.seed(12345)

def sys_from_parts(base_class, parts=None, gasparts=None, parameters=None,converter=None, extra=dict()):
  sys=base_class(convert_nbody=converter, **extra)
  sys.initialize_code()
  for param,value in parameters:
    err=sys.parameters.__setattr__(param,value)
  sys.commit_parameters()
  if parts is not None:
    sys.particles.add_particles(parts)
  if gasparts is not None:
    sys.gas_particles.add_particles(gasparts)
  if hasattr(sys,"start_viewer"):
    sys.start_viewer()
  return sys


Nstar=1000
total_mass=300 | units.MSun
Rscale=0.5 | units.parsec

conv = nbody_system.nbody_to_si(total_mass,Rscale)

parts=new_plummer_sphere(Nstar,convert_nbody=conv)
parts.radius=zero

fi=sys_from_parts(Fi,parts,converter=conv,parameters=[
                                   ("self_gravity_flag",False),
                                   ("adaptive_smoothing_flag",True),
                                   ("targetnn",7),
                                   ("nn_tol",0.01),
                                   ("epsilon_squared",zero)])

fi.commit_particles()

print fi.particles.radius[0:10]
