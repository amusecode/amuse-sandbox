"""
 run script for clustergas runs

"""

import numpy
import clustergas_gadget as clustergas
from amuse.support.units import units

numpy.random.seed(123491)

clustergas.clustergas(sfeff=0.3, 
                    Nstar=1000,
                    Ngas=100000,
                    Rscale=0.5 | units.parsec,
                    feedback_efficiency=0.01,
                    runid="test")


#clustergas.clustergas_restart(
#                   "test",
#                   250,newid='test_restart')

