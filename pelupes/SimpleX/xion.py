import numpy
import cPickle

from matplotlib import pyplot

from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.legacy.fi.interface import Fi
from amuse.support.data.core import Grid


from amuse.io import read_set_from_file

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

mu=1.| units.amu
muion =0.5 | units.amu
xtrans=0.06
mutrans=mu/(1+xtrans)
Tinit=100. | units.K
Ttrans=13500. | units.K 
Tion=13500. | units.K
rhoinit=0.001 | (units.amu / units.cm**3)
uinit=constants.kB * Tinit/mu
utrans=constants.kB * Ttrans/mutrans
uion=constants.kB * Tion/muion


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


def xion_from_u(u):
  xion=xtrans*(u-uinit)/(utrans-uinit)
  a=numpy.where(u>utrans)[0]
  xion[a]=xtrans+(1-xtrans)*(u[a]-utrans)/(uion-utrans)
  return xion

def plots(i):

  g=read_set_from_file('dump-%6.6i'%i,'amuse')

  r=((g.x**2+g.y**2+g.z**2)**0.5).value_in(units.kpc)
  v=((g.vx**2+g.vy**2+g.vz**2)**0.5).value_in(units.kms)
  cs=(g.u**0.5).value_in(units.kms)
#  xion=((g.u-uinit)/uion).value_in(units.none)
  xion=xion_from_u(g.u).value_in(units.none)
  rho=3./4/numpy.pi*8*g.mass/g.radius**3
  dens=(rho).value_in(units.amu/units.cm**3)
  pres=(g.u*rho).value_in(units.g/units.cm/units.s**2)
  mach=v/cs

  pyplot.figure(figsize=(6,6))
  pyplot.semilogy(r/15,xion,'r .')
  pyplot.semilogy(r/15,1-xion,'g .')
  pyplot.xlim((0.,1.))
  pyplot.ylim((1.e-6,1.))
  pyplot.savefig('xion-part-%6.6i.png'%i)

  aplot(i,'xion',((r/15,xion,'r'),(r/15,1-xion,'g')),
          xlim=(0.,1.),ylim=(1.e-6,1.))
  aplot(i,'pres',((r/15,pres,'r'),),
          xlim=(0.,1.),ylim=(1.e-17,1.e-14))
  aplot(i,'rho',((r/15,dens,'r'),),
          xlim=(0.,1.),ylim=(0.0001,0.01))
  aplot(i,'mach',((r/15,mach,'r'),),
          xlim=(0.,1.),ylim=(1.e-5,10.))

if __name__=="__main__":
#  plots(100)
  plots(300)
#  plots(500)

