import numpy
import cPickle

from matplotlib import pyplot,rc

from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.legacy.fi.interface import Fi
from amuse.io import read_set_from_file


from amuse.support.data import Grid
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

rc('font', **font)


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

def xion_from_u(u):
  xion=xtrans*(u-uinit)/(utrans-uinit)
  a=numpy.where(u>utrans)[0]
  xion[a]=xtrans+(1-xtrans)*(u[a]-utrans)/(uion-utrans)
  return xion

def plots(snaps):

  n=len(snaps)
  r=dict()
  xion=dict()

  for i in snaps:
    g=read_set_from_file('dump-%6.6i'%i,'amuse')
    r[i]=((g.x**2+g.y**2+g.z**2)**0.5).value_in(units.kpc)
#  xion=((g.u-uinit)/uion).value_in(units.none)
    xion[i]=xion_from_u(g.u).value_in(units.none)

  f=pyplot.figure(figsize=(8,8))
  for i,ii in enumerate(snaps):
    xa,ya=average(100,r[ii],xion[ii])
    subplot=f.add_subplot(n,1,i+1)

    subplot.semilogy(xa,ya,'r',linewidth=1.5)  
    subplot.semilogy(xa,1-ya,'g',linewidth=1.5)  
    subplot.set_ylabel('x, 1-x')  
    subplot.set_xlim( (0,15) )
    subplot.set_ylim( (1.e-4,1) )

    if i==n-1:
      subplot.set_xlabel('R (kpc)')  
    else:  
      subplot.set_xticklabels([])

    pyplot.figtext(0.75,0.85-i*0.83/n,'%3i Myr'%ii)
  pyplot.savefig('xion.eps')


if __name__=="__main__":
  plots([100,225,350,475])
