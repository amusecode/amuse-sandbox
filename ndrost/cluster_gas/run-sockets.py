"""
 run script for clustergas runs

"""

import numpy
import clustergas_gadget as clustergas
from amuse.units import units
from amuse.community.fi.interface import Fi
from amuse.community.octgrav.interface import Octgrav
from amuse.community.gadget2.interface import Gadget2
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPE
from SSEplus import SSEplus


#numpy.random.seed(123489)
#numpy.random.seed(123477)

clustergas.clustergas(sfeff=0.3, 
                    Nstar=100,
                    Ngas=10000,
                    Rscale=0.5 | units.parsec,
                    feedback_efficiency=0.1,
                    runid="demo-output",
		    dt_plot=0.01 | units.Myr,

                    grav_code=PhiGRAPE,
                    grav_code_extra=dict(mode='normal', channel_type='sockets'),

                    gas_code=Gadget2,
                    gas_code_extra=dict(mode='nogravity', number_of_workers=1,use_gl=False, channel_type='sockets'),

                    se_code=SSEplus,
                    se_code_extra=dict(channel_type='sockets'),

                    grav_couple_code=Fi,
                    grav_couple_code_extra=dict(channel_type='sockets')
)


#clustergas.clustergas_restart(
#                   "test",
#                   250,newid='test_restart')

