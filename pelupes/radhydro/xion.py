import numpy

from matplotlib import pyplot

from amuse.support.units import units
from amuse.support.units import constants

from amuse.ext.radial_profile import radial_profile

from amuse.support.io import read_set_from_file

def aplot(i, tag, xyc,xlim=None,ylim=None,ylabel=""):
  pyplot.figure(figsize=(8,8))
  for x,y,c in xyc:
    xa,ya=radial_profile(x,y,N=100)
    pyplot.semilogy(xa,ya,c)  
  pyplot.xlabel('L/Lbox')  
  pyplot.ylabel(ylabel)  
  if xlim is not None:
    pyplot.xlim(xlim)
  if ylim is not None:
    pyplot.ylim(ylim)
  pyplot.savefig(tag+'-%6.6i.png'%i)

def plots(i,label):

  g=read_set_from_file(label+'-%6.6i'%i,'amuse')

  gamma=5./3.
  mu=1.| units.amu

  r=((g.x**2+g.y**2+g.z**2)**0.5).value_in(units.kpc)
  v=((g.vx**2+g.vy**2+g.vz**2)**0.5).value_in(units.kms)
  cs=((gamma*(gamma-1.)*g.u)**0.5).value_in(units.kms)
  xion=g.xion
  T=((gamma-1)*g.u*mu/(1+xion)/constants.kB).value_in(units.K)
  rho=g.rho
  dens=(rho).value_in(units.amu/units.cm**3)
  pres=(g.u*rho).value_in(units.g/units.cm/units.s**2)
  mach=v/cs

  aplot(i,'xion',((r/15,xion,'r'),(r/15,1-xion,'g')),
          xlim=(0.,1.),ylim=(1.e-6,1.),ylabel="x, 1-x")
  aplot(i,'temp',((r/15,T,'r'),),
          xlim=(0.,1.),ylim=(0,1.e5),ylabel="T (K)")
  aplot(i,'pres',((r/15,pres,'r'),),
          xlim=(0.,1.),ylim=(1.e-17,1.e-14),ylabel="P (g/cm/s**2)")
  aplot(i,'rho',((r/15,dens,'r'),),
          xlim=(0.,1.),ylim=(0.0001,0.01),ylabel="density (amu/cm**3)")
  aplot(i,'mach',((r/15,mach,'r'),),
          xlim=(0.,1.),ylim=(1.e-5,10.),ylabel="mach number")
  aplot(i,'vel',((r/15,v,'r'),(r/15,cs,'g'),),
          xlim=(0.,1.),ylim=(0.01,100.),ylabel="bulk, sound speed (km/s)")

if __name__=="__main__":
  plots(10,'radhydro')

