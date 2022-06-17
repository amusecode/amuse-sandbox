import pickle

import numpy

from amuse.units import units

from amuse.community.sse.interface import SSE

from amuse.datamodel import Particles

from amuse.units.quantities import VectorQuantity

def stellar_remnant_state(star):
    return (10 <= star.stellar_type.value_in(units.stellar_type))* \
        (star.stellar_type.value_in(units.stellar_type) < 16)

def dimm_star(star):
    return (star.luminosity< 1.e-4 | units.LSun)

N=256

Mmin=0.1 | units.MSun
Mmax=100 | units.MSun

fac=numpy.log10(Mmax/Mmin)

mass=Mmin*10**(fac*numpy.arange(N)/(N-1))

parts=Particles(mass=mass)

se=SSE()

parts=se.particles.add_particles(parts)

data=[]

#~ done = False
#~ while not done:
  #~ se.evolve_model()  
  #~ data.append([parts.age, parts.time_step, parts.mass, parts.luminosity, parts.temperature])
  #~ done = stellar_remnant_state(parts)
  #~ print done
  #~ done = done.prod()
#~ 
#~ print len(data)
#~ 
for i,p in enumerate(se.particles):
  print(i,'/',N)
  pdata=[]
  done = False
  tlast=p.age
  fac=1.
  while not done:
    pdata.append([p.age, p.time_step, p.mass, p.luminosity, p.temperature])
    #~ done = stellar_remnant_state(p)
    done = dimm_star(p)
    timestep=p.time_step/4
    #~ if stellar_remnant_state(p):
      #~ timestep=timestep*fac
      #~ fac=min(fac*3,1024*2)
    #~ if timestep<1. | units.yr:
      #~ timestep=1.| units.yr
    p.evolve_for(timestep)
    tlast=p.age

  data.append(pdata)
#~ 


newdata=[]

for pdata in data:
  age=VectorQuantity.new_from_scalar_quantities(*[x[0] for x in pdata])
  timestep=VectorQuantity.new_from_scalar_quantities(*[x[1] for x in pdata])
  mass=VectorQuantity.new_from_scalar_quantities(*[x[2] for x in pdata])
  lum=VectorQuantity.new_from_scalar_quantities(*[x[3] for x in pdata])
  temp=VectorQuantity.new_from_scalar_quantities(*[x[4] for x in pdata])
  newdata.append([age,timestep,mass,lum,temp])


f=open("rundata.pkl","wb")
pickle.dump(newdata,f)
f.close()

for pdata in data:
  print(len(pdata))
  dt=[x[1] for x in pdata]
  dt=VectorQuantity.new_from_scalar_quantities(*dt)
  #~ print dt.min(),dt.max()

  lum=[x[3] for x in pdata]
  lum=VectorQuantity.new_from_scalar_quantities(*lum).in_(units.LSun)
  #~ print lum.min(),lum.max()



