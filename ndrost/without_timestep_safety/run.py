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
from amuse.community.ph4.interface import ph4
from SSEplus import SSEplus

#nice big star with 100 stars
#numpy.random.seed(123489)

#best guess inti's original seed
numpy.random.seed(123491)

clustergas.clustergas(sfeff=0.3, 
                    Nstar=10000,
                    Ngas=1000000,
                    Rscale=0.5 | units.parsec,
                    runid="sim-output",
                    feedback_efficiency=0.01,
		    dt_plot=0.05 | units.Myr,

		    #LGM
#                    grav_code=PhiGRAPE,
                    grav_code=ph4,
                    grav_code_extra=dict(mode='gpu', channel_type='ibis'),

		    #VU
                    gas_code=Gadget2,
                    gas_code_extra=dict(mode='nogravity', output_directory='output', number_of_workers=32, number_of_nodes=8,use_gl=False, channel_type='ibis', hostname='VU'),

                    #UVA
                    se_code=SSEplus,
                    se_code_extra=dict(channel_type='ibis'),

                    #DELFT
                    grav_couple_code=Octgrav,
                    grav_couple_code_extra=dict(channel_type='ibis')
)


#clustergas.clustergas_restart(
#                   "test",
#                   250,newid='test_restart')

