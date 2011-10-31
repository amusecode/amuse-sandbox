import numpy
import clustergas_gadget as clustergas
from amuse.support.units import units

#numpy.random.seed(123477)

clustergas.clustergas_restart(
                   "demo-output2",
                   160,
		   dt_plot=0.01 | units.Myr,
                   grav_code_extra=dict(mode='gpu', channel_type='sockets', hostname='localhost'),
                   gas_code_extra=dict(number_of_workers=8,number_of_nodes=1,use_gl=False, channel_type='ibis', hostname='local'),
                   se_code_extra=dict(channel_type='sockets', hostname='localhost'),
                   grav_couple_code_extra=dict(channel_type='sockets', hostname='localhost')

)
