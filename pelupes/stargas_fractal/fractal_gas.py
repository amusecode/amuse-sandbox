from amuse.test import amusetest

from amuse.ic.fractalcluster import new_fractal_cluster_model

from amuse.units import nbody_system
from amuse.units import units

from matplotlib import pyplot

N=100
Ng=10000

for fd in [1.6,1.8,2.,2.2,2.4,2.6,2.8,3.0]:

  stars = new_fractal_cluster_model(N=N,fractal_dimension=fd,random_seed=12345321)
  gas = new_fractal_cluster_model(N=Ng,fractal_dimension=fd,random_seed=12345321)
  
  f=pyplot.figure(figsize=(8,8))
  pyplot.plot(gas.x.number,gas.y.number,'g.')
  pyplot.plot(stars.x.number,stars.y.number,'r.')
  pyplot.savefig("stargas_fractal-%3.1f.png"%fd)
  
