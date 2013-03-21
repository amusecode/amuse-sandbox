import cPickle

from matplotlib import pyplot

from amuse.units.quantities import VectorQuantity

import numpy

new_vector=VectorQuantity.new_from_scalar_quantities

f=open('data','rb')
data=cPickle.load(f)
f.close()

t=new_vector(*data['time']).number

ek=new_vector(*data['kinetic_energy'])
ep=new_vector(*data['potential_energy'])
etot=ek+ep

err=abs((etot-etot[0])/etot[0])

mf=data['mf']
lagrangian_radii=dict()
for i,mf in enumerate(data['mf']):
  lagrangian_radii[mf]=map(lambda x: x[i].number, data['lagrangian_radii'])

rcore=map(lambda x: x.radius.number,data['core'])
rhocore=map(lambda x: x.density.number,data['core'])

wc=data['wallclock']

px=new_vector(*map(lambda x: x[0],data['total_momentum'])).number
py=new_vector(*map(lambda x: x[1],data['total_momentum'])).number
pz=new_vector(*map(lambda x: x[2],data['total_momentum'])).number

lx=new_vector(*map(lambda x: x[0],data['total_angular_momentum'])).number
ly=new_vector(*map(lambda x: x[1],data['total_angular_momentum'])).number
lz=new_vector(*map(lambda x: x[2],data['total_angular_momentum'])).number



pyplot.figure(figsize=(8,6))
pyplot.semilogy(t,rhocore,'b',lw=2)  
pyplot.xlabel('time')
pyplot.ylabel('core density')
#pyplot.show()
pyplot.savefig('rhocore.eps')


pyplot.figure(figsize=(8,6))
pyplot.semilogy(t,abs(lx-lx[0]),'r')
pyplot.semilogy(t,abs(ly-ly[0]),'r')
pyplot.semilogy(t,abs(lz-lz[0]),'r')
pyplot.xlabel('time')
pyplot.ylabel('L-L0')
#pyplot.ylim(-1.e-15,1.e-15)
#pyplot.show()
pyplot.savefig('l.eps')

pyplot.figure(figsize=(8,6))
pyplot.semilogy(t,abs(px),'r')
pyplot.semilogy(t,abs(py),'r')
pyplot.semilogy(t,abs(pz),'r')
pyplot.xlabel('time')
pyplot.ylabel('P')
#pyplot.ylim(-1.e-15,1.e-15)
pyplot.savefig('p.eps')

pyplot.figure(figsize=(8,6))
pyplot.semilogy(t,wc,'r')
pyplot.xlabel('simulation time')
pyplot.ylabel('wallclock time')
pyplot.savefig("wallclock.eps")

pyplot.figure(figsize=(8,6))
for mf,lr in lagrangian_radii.items():
  pyplot.semilogy(t,lr,'r')
pyplot.semilogy(t,rcore,'b',lw=2)  
pyplot.xlabel('time')
pyplot.ylabel('lagrangian, core radii')
pyplot.savefig('lr.eps')

pyplot.figure(figsize=(8,6))
pyplot.semilogy(t,err,'r')
pyplot.xlabel('time')
pyplot.ylabel('energy error')
pyplot.savefig('err.eps')


