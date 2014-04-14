from amuse.community.sse.interface import SSE

from amuse.units import units

from amuse.datamodel import Particles,Particle

def agb():

  p=Particles(1,mass=8. | units.MSun)
  
  sse=SSE()
  
  p=sse.particles.add_particles(p)[0]
  
  time=[] | units.Myr
  mass=[] | units.MSun
  
  time.append(p.age)
  mass.append(p.mass)
  
  print time
  
  while p.mass > (1.438 | units.MSun): 
    p.evolve_for( 1.e3 | units.yr)
    time.append(p.age)
    mass.append(p.mass)
  
  print time[-1]
    
  from matplotlib import pyplot  
  
  pyplot.plot(time.value_in(units.Myr),mass.value_in(units.MSun))
  pyplot.show()
  
def agb_test():
  
  p=Particle(mass=8. | units.MSun)
  sse=SSE()
  
  p=sse.particles.add_particle(p)
  sse.evolve_model(41.9 | units.Myr)
  
  p=Particle(mass=8. | units.MSun,age=p.age)
  sse=SSE()
  p=sse.particles.add_particle(p)
    
  time=[] | units.Myr
  mass=[] | units.MSun
  
  time.append(p.age)
  mass.append(p.mass)
  
  print time
  
  while p.mass > (1.438 | units.MSun): 
    p.evolve_for( 1.e3 | units.yr)
    time.append(p.age)
    mass.append(p.mass)
  
  print time[-1]
  
  from matplotlib import pyplot  
  
  pyplot.plot(time.value_in(units.Myr),mass.value_in(units.MSun))
  pyplot.show()  

if __name__=="__main__":
  agb()
