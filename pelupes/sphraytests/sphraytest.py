import os.path
import numpy

from amuse.community.sphray.interface import SPHRay
from amuse.units import units,constants
from amuse.datamodel import Particles,create_particle_set

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

def aplot(i, tag, xyc,xlim=None,ylim=None,ylab=''):
  pyplot.figure(figsize=(6,6))
  for x,y,c in xyc:
    xa,ya=average(100,x,y)
    pyplot.semilogy(xa,ya,c)
  pyplot.xlabel('L/Lbox')
  pyplot.ylabel(ylab)
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
    return create_particle_set(mass=mass | (10**10*units.MSun), h_smooth=hsml | (units.kpc), 
        x=x | (units.kpc), y=y| (units.kpc), z=z| (units.kpc), rho=rho | ((10**10*units.MSun) /(units.kpc)**3),
        xion=xe, u=u| (10**5 * units.cm/units.s)**2)
    
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
    return create_particle_set( luminosity=1.e50*numpy.array(L) | ( units.s**-1), 
        x=numpy.array(x) | (units.kpc), y=numpy.array(y)| (units.kpc), 
              z=numpy.array(z)| (units.kpc), SpcType=numpy.array(spctype))


if __name__=="__main__":
  rad = SPHRay(redirection="none")#,channel_type="ibis", hostname="paddegat")
  
  rad.parameters.isothermal_flag=False
  rad.parameters.number_of_rays=10000

  print rad.parameters

  input_file = 'sphray_glass_L13.2_N64'
  gas=read_gas_file(input_file)
#  gas.u=gas.u*100

  input_file = os.path.join(rad.data_directory(), 'test2_sources_001.1')
  src = read_src_file(input_file)

  rad.src_particles.add_particles(src)
  rad.gas_particles.add_particles(gas)
   
  t1=time.time()
  rad.evolve_model(15. | units.Myr)
  rad.evolve_model(30. | units.Myr)
  t2=time.time()
 
  x=rad.gas_particles.x.value_in(units.kpc)
  y=rad.gas_particles.y.value_in(units.kpc)
  z=rad.gas_particles.z.value_in(units.kpc)
  xs=src.x.value_in(units.kpc)
  ys=src.y.value_in(units.kpc)
  zs=src.z.value_in(units.kpc)
  r=((x-xs[0])**2+(y-ys[0])**2+(z-zs[0])**2)**0.5

  xion=rad.gas_particles.xion
  
  gamma=1.66667
  mu=1.| units.amu
  t=(gamma-1)*mu/(1.+xion)/constants.kB*rad.gas_particles.u
  t=t.value_in(units.K)

  aplot(0,'xion',((r/6.6,xion,'r'),(r/6.6,1-xion,'g')),
            xlim=(0.,1.),ylim=(1.e-6,1.),ylab='x, 1-x')
  aplot(0,'temp',((r/6.6,t,'r'),),
            xlim=(0.,1.),ylim=(10,1.e5),ylab='T (K)')
  
  print t2-t1
  
    
  
