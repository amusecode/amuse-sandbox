import os.path
import numpy
from amuse.test.amusetest import TestWithMPI

from amuse.community.sphray.interface import SPHRayInterface
from amuse.units import units
from amuse.datamodel import Particles

from amuse.io import read_set_from_file

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

import time

def average(N,r,dat):
  n=len(r)
  a=numpy.argsort(r)
  i=0
  r_a=[]
  dat_a=[]
  while i < n:
    ra=r[a[i:i+N]].sum()/min(n-i,N)
    da=dat[a[i:i+N]].sum()/min(n-i,N)
    r_a.append(ra)
    dat_a.append(da)
    i=i+N
  return numpy.array(r_a),numpy.array(dat_a)

def aplot(i, tag, xyc,xlim=None,ylim=None):
  pyplot.figure(figsize=(6,6))
  for x,y,c in xyc:
    xa,ya=average(100,x,y)
    pyplot.semilogy(xa,ya,c)
  pyplot.xlabel('L/Lbox')
  pyplot.ylabel('x, 1-x')
  if xlim is not None:
    pyplot.xlim(xlim)
  if ylim is not None:
    pyplot.ylim(ylim)
  pyplot.savefig(tag+'-%6.6i.png'%i)



def read_gas_file(filename):
    p=read_set_from_file(filename,'amuse')    
    mass=p.mass.number
    hsml=p.smoothing_length.number
    x=p.x.number
    y=p.y.number
    z=p.z.number
    rho=p.rho.number
    u=p.internal_energy.number
    xe=numpy.zeros_like(x)
    return mass, hsml, x, y, z, rho, xe, u

def read_src_file(filename):
    f=open(filename)
    lines=f.readlines()
    f.close()
    L=[]
    x=[]
    y=[]
    z=[]
    spctype=[]
    for line in lines:
      l=line.split()
      if len(l) == 9:
        L.append(float(l[6]))
        x.append(float(l[0]))
        y.append(float(l[1]))
        z.append(float(l[2]))
        spctype.append(float(l[7]))
    return numpy.array(L), numpy.array(x), numpy.array(y), numpy.array(z), numpy.array(spctype)


instance = SPHRayInterface(channel_type="ibis", hostname="paddegat")
instance.initialize_code()

instance.set_data_directory(instance.data_directory())
instance.set_output_directory(instance.output_directory())  

instance.set_isothermal(1)

instance.commit_parameters()

input_file = 'sphray_glass_L13.2_N64'
mass, hsml, x, y, z, rho, xe, u = read_gas_file(input_file)

u=u*100.

gas, errors = instance.new_gas_particle(mass, hsml, x, y, z, rho, xe, u)

input_file = os.path.join(instance.data_directory(), 'test2_sources_001.1')
L, xs, ys, zs, spctype = read_src_file(input_file)

src, errors = instance.new_src_particle(L, xs, ys, zs, spctype)

print len(gas),len(src)

instance.commit_particles()

t1=time.time()
instance.evolve_model(30./978)
t2=time.time()

mass, hsml, x, y, z, rho, xe, u , error = instance.get_state_gas(gas)

r=((x-xs[0])**2+(y-ys[0])**2+(z-zs[0])**2)**0.5


#pyplot.figure(figsize=(8,6))
#pyplot.plot(r,xe,'r.')
#pyplot.savefig('r-xe.png')
#pyplot.figure(figsize=(8,6))
#pyplot.plot(r,rho,'r.')
#pyplot.savefig('r-rho.png')
#pyplot.figure(figsize=(8,6))
#pyplot.plot(r,u,'r.')
#pyplot.savefig('r-u.png')

aplot(0,'xion',((r/6.6,xe,'r'),(r/6.6,1-xe,'g')),
          xlim=(0.,1.),ylim=(1.e-6,1.))
aplot(0,'xion',((r/6.6,u,'r')),
          xlim=(0.,1.),ylim=(1.e-6,1.))

print t2-t1
