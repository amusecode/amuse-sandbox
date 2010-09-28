# estimate the mean maximum mass star for a cluster with a Salpeter IMF
# this illustrates that the units won't get in the way of calculations
# other than that the units don't play a role.

#import numpy and set random seed
import numpy
numpy.random.seed(123456)

#import AMUSE units and Salpeter IMF
from amuse.support.units import units
from amuse.ext.salpeter import SalpeterIMF

# calculate the maximum stellar mass for Nset realizations of 
# an N star cluster
N=1000
Nset=100

maxmasses=[]|units.MSun
for i in range(Nset):
  tm,m=SalpeterIMF(mass_max = 100. | units.MSun).next_set(N)
  maxmasses.append(max(m))
  
print "median:",maxmasses.median()
print "mean:", maxmasses.mean()
# result:
#  median: 19.6164811985 MSun
#  mean: 26.1161949491 MSun

print "stddev:", numpy.std(maxmasses)
#result:
# stddev: 19.2701229918 1.98892e+30 * kg

# this can be fixed
print "stddev:", numpy.std(maxmasses).in_(units.MSun)
# result:
# stddev: 19.2701229918 MSun

# also:
print "median:",numpy.median( maxmasses).in_(units.MSun)
# result:
# median: 19.6164811985 MSun

# however:
print "mean:", numpy.mean( maxmasses)
# gives an error!



