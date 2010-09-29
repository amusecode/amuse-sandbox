# estimate the mean maximum mass star for a cluster with a Salpeter IMF
# this illustrates that the units won't get in the way of calculations
# other than that the units don't play a role.

#import numpy and set random seed
import numpy
numpy.random.seed(123456)

#import AMUSE units and Salpeter IMF
from amuse.support.units import units
from amuse.ext.salpeter import new_salpeter_mass_distribution

# calculate the maximum stellar mass for Nset realizations of 
# an N star cluster
N = 1000
Nset = 100

maxmasses = []|units.MSun
for i in range(Nset):
    maxmasses.append(max(new_salpeter_mass_distribution(N, mass_max = 100. | units.MSun)))

print "    using AMUSE VectorQuantity functions"
print "mean:  ", maxmasses.mean()
print "stddev:", maxmasses.std()
print "median:", maxmasses.median()
print
print "    using numpy functions"
print "mean:  ", numpy.mean( maxmasses)
print "stddev:", numpy.std(maxmasses)
print "median:", numpy.median(maxmasses)

# should result in something like:

# mean:   27.4915750164 MSun
# stddev: 19.7149800906 MSun
# median: 21.0983403429 MSun
# mean:   27.4915750164 MSun
# stddev: 19.7149800906 MSun
# median: 21.0983403429 1.98892e+30 * kg

# the conversion of the unit of the final output can be fixed:
print "median:", numpy.median(maxmasses).in_(units.MSun)
# median: 21.0983403429 MSun
