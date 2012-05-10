import os.path
import numpy as np
from amuse.test.amusetest import TestWithMPI

from amuse.community.sphray.interface import SPHRayInterface
from amuse.units import units
from amuse.datamodel import Particles

from amuse.io import read_set_from_file

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

import time

class UserDefinedType:
  pass







# average quantity with an equal number of particles in each bin
#-----------------------------------------------------------------
def average(NinBin, r, dat):

  bins = UserDefinedType()
  bins.r = []
  bins.d = []

  N = len(r)
  indx = np.argsort(r)
  i = 0

  while i < N:
    r_av = np.sum(   r[indx[i:i+NinBin]] ) / min(N-i,NinBin)
    d_av = np.sum( dat[indx[i:i+NinBin]] ) / min(N-i,NinBin)
    bins.r.append(r_av)
    bins.d.append(d_av)
    i = i + NinBin

  bins.r = np.array( bins.r )
  bins.d = np.array( bins.d )

  return bins.r, bins.d


# wrapper to matplotlib plot
#-----------------------------------------------------------------
def aplot(i, tag, xyc, xlim=None, ylim=None):

  pyplot.figure(figsize=(6,6))

  for x,y,c in xyc:
    xa,ya = average(1000, x, y)
    pyplot.semilogy(xa, ya, c)

  pyplot.xlabel('L/Lbox')
  pyplot.ylabel('x, 1-x')

  if xlim is not None:
    pyplot.xlim(xlim)

  if ylim is not None:
    pyplot.ylim(ylim)

  pyplot.savefig(tag+'-%6.6i.png'%i)


# Amuse I/O density field
#-----------------------------------------------------------------
def read_gas_file(filename):
    p = read_set_from_file(filename,'amuse')    
    mass=p.mass.number
    hsml=p.smoothing_length.number
    x=p.x.number
    y=p.y.number
    z=p.z.number
    rho=p.rho.number
    u=p.internal_energy.number
    xe=np.zeros_like(x)
    return mass, hsml, x, y, z, rho, xe, u


# Amuse I/O sources
#-----------------------------------------------------------------
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
    return np.array(L), np.array(x), np.array(y), np.array(z), np.array(spctype)



# Begin Test Code
#========================================================================


# initialize instance
#------------------------------------------------------------
#instance = SPHRayInterface(redirection='none')
instance = SPHRayInterface()
instance.initialize_code()


# commit instance parameters
#------------------------------------------------------------
instance.set_data_directory(instance.data_directory())
instance.set_output_directory(instance.output_directory())  

instance.set_isothermal(1)

instance.commit_parameters()


# read and import gas and sources
#------------------------------------------------------------
input_file = 'sphray_glass_L13.2_N64'
mass, hsml, x, y, z, rho, xe, u = read_gas_file(input_file)

u=u*100.

gas, errors = instance.new_gas_particle(mass, hsml, x, y, z, rho, xe, u)

input_file = os.path.join(instance.data_directory(), 'test1_sources_001.1')
L, xs, ys, zs, spctype = read_src_file(input_file)

src, errors = instance.new_src_particle(L, xs, ys, zs, spctype)

print 'ngas: ', len(gas)
print 'nsrc: ', len(src)

instance.commit_particles()

t1=time.time()
instance.evolve_model(30./978)
t2=time.time()
print 'time to evolve model: ', t2-t1


# report results
#------------------------------------------------------------

mass, hsml, x, y, z, rho, xe, u , error = instance.get_state_gas(gas)

print
print 'min/max: x:    ', x.min(), x.max()
print 'min/max: y:    ', y.min(), y.max()
print 'min/max: z:    ', z.min(), z.max()
print 'min/max: mass: ', mass.min(), mass.max()
print 'min/max: hsml: ', hsml.min(), hsml.max()
print 'min/max: rho:  ', rho.min(), rho.max()
print 'min/max: xe:   ', xe.min(), xe.max()
print 'min/max: u:    ', u.min(), u.max()
print



r = ( (x-xs[0])**2 + \
      (y-ys[0])**2 + \
      (z-zs[0])**2 )**0.5


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



