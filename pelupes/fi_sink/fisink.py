from amuse.units import nbody_system
from amuse.units import units

import amuse.datamodel as core
from amuse.community.fi.interface import Fi

from amuse.ext.sink import SinkParticles

"""
SinkFi

fi with sink particles

initialize with:
fi=Fi(density_threshold=maximum_density)

this will make a Fi SPH code where during evolve particles which exceed the maximum density are either
converted to sink particles or accreted by existing sink particles. Low density particles are NOT
accreted.

- writing out and restarting:

write out gas and dm sets as usual, in addition:
write_set_to_file(fi.sink,filename, file_format)

- when restarting:

read in all particle sets, and additionally the sink set, and initialize:

fi=Fi(density_threshold=maximum_density)
<set parameters>
fi.gas_particles.add_particles(gas)
fi.dm_particles.add_particles(dm)
fi.reinit_sink(sink)


"""



class SinkFi(Fi):
    def __init__(self, *args, **kargs):
        Fi.__init__(self, *args, **kargs)
        self.sink=None
        if not kargs.has_key('density_threshold'):
          raise Exception("provide density threshold")
        self.density_threshold=kargs['density_threshold']
    
    def evolve_model(self, *args, **kargs):
        self.parameters.stopping_condition_maximum_density = self.density_threshold
        density_limit_detection = self.stopping_conditions.density_limit_detection
        density_limit_detection.enable()
        self.overridden().evolve_model(*args,**kargs)        
        while density_limit_detection.is_set():
          highdens = density_limit_detection.particles().copy_to_memory()
          self.gas_particles.remove_particles(highdens)
          if self.sink is not None:
            self.sink.accrete(highdens)
            highdens_in_code = self.dm_particles.add_particles(highdens)
            self.sinks.add_sinks(highdens_in_code)
          else:
            highdens_in_code = self.dm_particles.add_particles(highdens)
            self.sink=SinkParticles(highdens_in_code)
          self.overridden().evolve_model(*args,**kargs)
          
    def reinit_sink(self,sinks):
        self.sink=SinkParticles(sinks.get_intersecting_subset_in(self.dm_particles),sink_radius=sinks.sink_radius)


