
from amuse.units import units,nbody_system

from amuse.io import read_set_from_file

from amuse.community.huayno.interface import Huayno

from amuse.datamodel import Particles, ParticlesOverlay

from matplotlib import pyplot

def read_encounters(filename):
  encounters = read_set_from_file(filename, 'txt', 
    attribute_names=['t_enc', 'mass', 'vx_enc','vy_enc', 'vz_enc', 'x_start', 'y_start', 'z_start', 'timestep' ],
    attribute_types=[units.yr, units.MSun,  units.AU/units.yr, units.AU/units.yr, units.AU/units.yr, units.AU, units.AU, units.AU, units.yr],
    header_prefix_string='"' )
  return encounters  

class EncounterSystem(Huayno):
    def __init__(self,converter,filename='generate_encounter_1Gyr_simple.txt'):
        Huayno.__init__(self,converter)
        encounters=read_encounters(filename)
        encounters.t_start=encounters.t_enc-3*encounters.timestep
        encounters.t_end=encounters.t_enc+3*encounters.timestep
        self.encounters=encounters.sorted_by_attribute("t_start")
        self.particles=ParticlesOverlay(self.overridden().particles)
        self.current_index=0
        
    def next_finished_encounter(self,tend):
        if len(self.particles)==0:
          return tend,None
        part=self.particles.sorted_by_attribute("t_end")
        if part[0].t_end>tend:
          return tend,None
        else:
          return part[0].t_end,part[0]  
        
    def next_new_encounter(self,tend):    
        if self.current_index==len(self.encounters):
          return tend,None
        if self.encounters[self.current_index].t_start>tend:
          return tend,None
        else:
          part=self.encounters[self.current_index]
          part.x=part.x_start
          part.y=part.y_start
          part.z=part.z_start
          part.vx=part.vx_enc
          part.vy=part.vy_enc
          part.vz=part.vz_enc
          self.current_index+=1
          return part.t_start,part
        
    def evolve_model(self,tend):
        tnext=self.model_time
        while tnext<tend:
          tnew,new=self.next_new_encounter(tend)
          tfinish,finish=self.next_finished_encounter(tend)
          tnext=min(tend,tnew,tfinish)
          print len(self.particles),self.model_time.in_(units.Myr),tend,tnew,tfinish
          self.overridden().evolve_model(tnext)          
          if tnew==tnext and new is not None:
            print "add"
            print new.t_start
            self.particles.add_particle(new)
          if tfinish==tnext and finish is not None:  
            print "remove"
            self.particles.remove_particle(finish)

if __name__=="__main__":

  conv=nbody_system.nbody_to_si(1000 | units.AU,1. | units.kms)
  
  enc=EncounterSystem(conv)

  enc.begin_time=enc.encounters[0].t_start

  f=pyplot.figure(figsize=(8,8))
#  pyplot.ion()
#  pyplot.show()

  tend=1| units.Gyr
  dt=0.001| units.Myr

  tnow=0 | units.Myr

  i=0
  while tnow<tend:
    tnow+=dt
    enc.evolve_model(tnow+dt)
    print tnow,enc.model_time

    f=pyplot.figure(figsize=(8,8))
    pyplot.plot( 0.,0.,'b+')
    if len(enc.particles)>0:
      pyplot.plot(enc.particles.x.value_in(units.AU),enc.particles.position.lengths().value_in(units.AU),'ro')
    pyplot.xlim(-1.e6,1.e6)
    pyplot.ylim(0,2.e6)
    pyplot.xlabel('x [AU]')
    pyplot.ylabel('r [AU]')

    pyplot.savefig("frames/frame-%6.6i.png"%i)
    i+=1
    pyplot.clf()
    pyplot.close(f)
