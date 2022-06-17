import numpy
from PIL import Image

from amuse.units import units

import pickle

from amuse.units.quantities import VectorQuantity

from matplotlib import pyplot


def convert_from_byte(x,mn,mx):
  return mn*10**(x/255.*numpy.log10(mx/mn))

Nplot=32
Ndata=32

f=open("minmax%i.pkl"%Ndata,"rb")
dtmin,dtmax,Mmin,Mmax,Lmin,Lmax,Tmin,Tmax=pickle.load(f, encoding="latin1")
f.close()

print(dtmin,dtmax)
print(Mmin,Mmax)
print(Lmin,Lmax)
print(Tmin,Tmax)


f=open("rundata.pkl","rb")
data=pickle.load(f, encoding="latin1")
f.close()

for pdata in data[::len(data)//Nplot]:
  age=pdata[0]
  timestep=pdata[1]
  mass=pdata[2]
  lum=pdata[3]
  temp=pdata[4]
  
  x=age.value_in(units.Myr)
  x[0]+=0.01
  y=lum.value_in(units.LSun)
  #~ y=timestep.value_in(units.Myr)
  #~ y=temp.value_in(units.K)
  #~ y=mass.value_in(units.MSun)
    
  pyplot.loglog(x, y,'+')


f=open("data%i.pkl"%Ndata,"rb")
data=pickle.load(f, encoding="latin1")
f.close()

for pdata in data[::len(data)//Nplot]:
  age=pdata[0]
  timestep=pdata[1]
  mass=pdata[2]
  lum=pdata[3]
  temp=pdata[4]

  x=age.value_in(units.Myr)
  x[0]+=0.01
  y=lum.value_in(units.LSun)
  #~ y=timestep.value_in(units.Myr)
  #~ y=temp.value_in(units.K)
    
  #~ pyplot.loglog(x, y,'+')
  pyplot.loglog(x,y)

data=Image.open("data%i.png"%Ndata)

rdata=[]

for mi in range(Ndata):
  age=numpy.zeros(Ndata) | units.Myr
  mass=numpy.zeros(Ndata) | units.MSun
  lum=numpy.zeros(Ndata) | units.LSun
  temp=numpy.zeros(Ndata) | units.K

  tnew=0| units.Myr
  for ti in range(Ndata):
    age[ti]=tnew
    l,t,dt,m=data.getpixel((mi,ti))
    dt=convert_from_byte(dt,dtmin,dtmax)
    m=convert_from_byte(m,Mmin,Mmax)
    l=convert_from_byte(l,Lmin,Lmax)
    t=convert_from_byte(t,Tmin,Tmax)
    mass[ti]=m
    lum[ti]=l
    temp[ti]=t
    tnew=tnew+dt

  rdata.append([age,mass,lum,temp])

for pdata in rdata[::Ndata//Nplot]:
  age=pdata[0]
  mass=pdata[1]
  lum=pdata[2]
  temp=pdata[3]

  x=age.value_in(units.Myr)
  x[0]+=0.01
  y=lum.value_in(units.LSun)
  #~ y=timestep.value_in(units.Myr)
  #~ y=temp.value_in(units.K)
    
  pyplot.loglog(x, y,'p')
  #~ pyplot.loglog(x,y)

pyplot.show()
