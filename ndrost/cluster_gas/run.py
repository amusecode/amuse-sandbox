"""
 run script for clustergas runs

"""

import numpy
import clustergas_gadget as clustergas
from amuse.units import units

numpy.random.seed(123489)

clustergas.clustergas(sfeff=0.3, 
                    Nstar=100,
                    Ngas=10000,
                    Rscale=0.5 | units.parsec,
                    feedback_efficiency=0.01,
                    runid="demo-output")


#clustergas.clustergas_restart(
#                   "test",
#                   250,newid='test_restart')

