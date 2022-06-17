import numpy
from PIL import Image

from amuse.units import units

import pickle

from amuse.units.quantities import VectorQuantity

from matplotlib import pyplot

from resample import resample_tracking_extrema, resample_minimize

f=open("rundata.pkl","rb")
data=pickle.load(f)
f.close()

Ntarget=32
Nmass=32

min_timestep=1.e9 | units.Myr
max_timestep=0. | units.Myr
min_mass=1.e9 | units.MSun
max_mass=0. | units.MSun
min_lum=1.e9 | units.LSun
max_lum=0. | units.LSun
min_temp=1.e9 | units.K
max_temp=0. | units.K

dataeq=[]

for pdata in data[::len(data)//Nmass]:
  age=pdata[0]
  timestep=pdata[1]
  mass=pdata[2]
  lum=pdata[3]
  temp=pdata[4]
  N=len(age)
  
  
  dn=N/(1.*Ntarget)
  #~ print dn
  
  print("N:",N)
  
  ind=numpy.arange(N)

  #~ ind=resample_tracking_extrema(Ntarget,ind, numpy.log10(lum.value_in(units.LSun)))
  age_=age.copy()
  age_[0]=0.01 | units.Myr

  ind=resample_minimize(Ntarget,ind, numpy.log10(age_.value_in(units.Myr)),
                 numpy.log10(lum.value_in(units.LSun)),
                 numpy.log10(temp.value_in(units.K)),
                 numpy.log10(mass.value_in(units.MSun)))


  if N>Ntarget:
    assert len(ind)==Ntarget
  else:
    assert len(ind)==N
  
  x=age.value_in(units.Myr)
  x[0]+=0.01
  y=lum.value_in(units.LSun)
  #~ y=timestep.value_in(units.Myr)
  y2=temp.value_in(units.K)
  #~ y=mass.value_in(units.MSun)
    
  pyplot.loglog(x[ind],y[ind])
  pyplot.loglog(x[ind],y[ind],'p')
  pyplot.loglog(x, y,'+')
#~ 
  #~ pyplot.loglog(x[ind],y2[ind])
  #~ pyplot.loglog(x[ind],y2[ind],'p')
  #~ pyplot.loglog(x, y2,'+')



  age=age[ind]
  timestep=age[1:]-age[:-1]
  mass=mass[ind]
  lum=lum[ind]
  temp=temp[ind]
  
  
  #~ print timestep.min()

  min_timestep=min(timestep.min(), min_timestep)
  max_timestep=max(timestep.max(), max_timestep)

  min_mass=min(mass.min(), min_mass)
  max_mass=max(mass.max(), max_mass)

  min_lum=min(lum.min(), min_lum)
  max_lum=max(lum.max(), max_lum)

  min_temp=min(temp.min(), min_temp)
  max_temp=max(temp.max(), max_temp)

  dataeq.append([age,timestep,mass,lum,temp])

f=open("data%i.pkl"%Nmass,"wb")
pickle.dump(dataeq,f)
f.close()

min_mass=0.1 | units.MSun # since then initial masses are exact
max_mass=100. | units.MSun # since then initial masses are exact
min_lum=1.e-4 | units.LSun # floor in L
min_temp=2.e3 | units.K
max_temp=2.2e6 | units.K

print(min_timestep,max_timestep)
print(min_mass,max_mass)
print(min_lum,max_lum)
print(min_temp,max_temp)

f=open("minmax%i.pkl"%Nmass,"wb")
pickle.dump([min_timestep,max_timestep,
              min_mass,max_mass,
              min_lum,max_lum,
              min_temp,max_temp],f)
f.close()

mass_byt=numpy.zeros((Ntarget,Nmass),'uint8')
lum_byt=numpy.zeros((Ntarget,Nmass),'uint8')
temp_byt=numpy.zeros((Ntarget,Nmass),'uint8')
timestep_byt=numpy.zeros((Ntarget,Nmass),'uint8')


for i,pdata in enumerate(dataeq):
  age, timestep, mass, lum, temp=pdata
  
  mass=numpy.clip(mass.value_in(units.MSun),min_mass.value_in(units.MSun),max_mass.value_in(units.MSun)) | units.MSun
  mass=255*numpy.log10(mass/min_mass)/numpy.log10(max_mass/min_mass)
  mass=numpy.array(numpy.round(mass),'uint8')
  
  lum=numpy.clip(lum.value_in(units.LSun), 
                  min_lum.value_in(units.LSun),max_lum.value_in(units.LSun)) | units.LSun
  lum=255*numpy.log10(lum/min_lum)/numpy.log10(max_lum/min_lum)
  lum=numpy.array(numpy.round(lum),'uint8')

  temp=numpy.clip(temp.value_in(units.K), 
                  min_temp.value_in(units.K),max_temp.value_in(units.K)) | units.K
  temp=255*numpy.log10(temp/min_temp)/numpy.log10(max_temp/min_temp)
  temp=numpy.array(numpy.round(temp),'uint8')

  mass_byt[:len(mass),i]=mass
  lum_byt[:len(lum),i]=lum
  temp_byt[:len(temp),i]=temp

  tnow=age[0]
  timestep=numpy.zeros(len(age),'uint8')
  for k,ttarget in enumerate(age[1:]):
    dt=ttarget-tnow
    dt=numpy.clip(dt,min_timestep,max_timestep)
    dti=255*numpy.log10(dt/min_timestep)/numpy.log10(max_timestep/min_timestep)
    timestep[k]=numpy.int(numpy.floor(dti))
    tnow=tnow+min_timestep*10**(timestep[k]/255.*numpy.log10(max_timestep/min_timestep))
    
  timestep_byt[:len(timestep),i]=timestep

a=Image.fromarray(mass_byt)
r=Image.fromarray(lum_byt)
g=Image.fromarray(temp_byt)
b=Image.fromarray(timestep_byt)

argb=Image.merge("RGBA",(r,g,b,a))

argb.save("data%i.png"%Nmass)

pyplot.show()
    
pyplot.imshow(numpy.transpose(lum_byt),origin='lower')
pyplot.show()
