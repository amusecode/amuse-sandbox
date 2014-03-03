from amuse.units import nbody_system
from amuse.units import units

from amuse.datamodel import Particle,Particles
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


note hard coded factor 2 (appropiate for fi: sink radius is smoothing kernel size)

"""

class SinkFi(Fi):
    def __init__(self, *args, **kargs):
        self.sink_particles=Particles()
        if not kargs.has_key('density_threshold'):
          raise Exception("provide density_threshold")
        self.density_threshold=kargs.pop('density_threshold')
        if not kargs.has_key('merge_radius'):
          raise Exception("provide merge_radius")
        self.merge_radius=kargs.pop('merge_radius')
        self.verbose=kargs.pop("verbose", False)
        Fi.__init__(self, *args, **kargs)

    def merge_sinks(self):
        if self.verbose: print "identify groups.." 
        ccs=self.sink_particles.copy().connected_components(threshold=self.merge_radius)
        if len(ccs):
            if self.verbose:  print "merging sink sets... "
        nmerge=0
        newsinks=Particles()
        for cc in ccs:
            if len(cc) >1:
              nmerge+=1
              if len(cc) > 3: 
                print "warning: large merge(than 3) ", len(cc)
              cc=cc.copy()
              self.sink_particles.remove_particles(cc)
              new=cc.sorted_by_attribute('mass')[-1].empty_copy()
              try:
                new.radius=cc.sink_radius.max()
              except Exception as ex:
                print cc.__class__.__name__
                print dir(cc)
                print cc.sink_radius.__class__.__name__
                print dir(cc.sink_radius)
                raise ex 
              new.mass=cc.total_mass()
              new.position=cc.center_of_mass()
              new.velocity=cc.center_of_mass_velocity()
              cc.move_to_center()
              new.angular_momentum=cc.total_angular_momentum()
              new.sink_radius=new.radius
              newsinks.add_particle(new)
        if len(newsinks)>0:
            new_in_code = self.dm_particles.add_particles(newsinks)
            self.sink_particles.add_sink(new_in_code)
        if nmerge>0:
          if self.verbose: print "nmerg found: ",nmerge,
        if self.verbose: print "...done"

    def commit_parameters(self):
        self.parameters.stopping_condition_maximum_density = self.density_threshold
        self.overridden().commit_parameters()

    def evolve_model(self, *args, **kargs):
        density_limit_detection = self.stopping_conditions.density_limit_detection
        density_limit_detection.enable()
        self.overridden().evolve_model(*args,**kargs)
        while density_limit_detection.is_set():
          if self.verbose: print "processing high dens particles..."
          highdens=self.gas_particles.select_array(lambda rho:rho> self.density_threshold,["rho"])
          candidate_sinks=highdens.copy()
          if len(candidate_sinks)==0:
            print "WARNING: unexpected no highdens"
            raise Exception("I better stop")
          candidate_sinks.radius*=2 
          if len(self.sink_particles)>0:
            if self.verbose: print "accreting and new sinks..."
            self.sink_particles.accrete(candidate_sinks)
            self.gas_particles.remove_particles(highdens)
            if len(candidate_sinks)>0:
              newsinks_in_code = self.dm_particles.add_particles(candidate_sinks)
              self.sink_particles.add_sinks(newsinks_in_code)
          else:
            if self.verbose: print "new sinks..."
            self.gas_particles.remove_particles(highdens)
            newsinks_in_code = self.dm_particles.add_particles(candidate_sinks)
            self.sink_particles=SinkParticles(newsinks_in_code,looping_over="sources")
          if self.verbose: print "..done"
          self.overridden().evolve_model(*args,**kargs)
        if len(self.sink_particles)>1: self.merge_sinks()        
  
    def reinit_sink(self,sink_particles):
        sinks_in_code=sink_particles.get_intersecting_subset_in(self.dm_particles)
        if len(sinks_in_code) != len(sink_particles):
          raise Exception("all sinks need to be in dm_particles")
        self.sink_particles=SinkParticles(sinks_in_code,sink_radius=sink_particles.sink_radius,looping_over="sources")
