import os
import sys
import numpy

from matplotlib import pyplot
    
from amuse.units import nbody_system
from amuse.units import units
from amuse.io import write_set_to_file    
    
from amuse.ext.radial_profile import radial_density

from amuse.io import read_set_from_file

fgas=0.5

def radial_plot(rad_dens,label="none"):
  c=['g','b','y','c','r:','g:']

  f=pyplot.figure(figsize = (8, 6))
  pyplot.xlabel(r'radius')
  pyplot.ylabel(r'density')
  for i in range(len(rad_dens)):
    pyplot.loglog(rad_dens[i][0],rad_dens[i][1],c[i%len(c)])
  
  ascl=1/1.695
  ra,dens=rad_dens[0]
  pyplot.loglog(ra, fgas* 3./4./numpy.pi/ascl**3/(1+(ra**2/ascl**2))**(5./2),'r')
  pyplot.xlim(0.1,10)
  pyplot.ylim(1.e-7,10)
  pyplot.savefig(label+'-rad_dens.eps')
  pyplot.close(f)

rad_dens_gas=[]
rad_dens_stars=[]

for i in range(2,11,2):
    print i
    gas=read_set_from_file('gas-%6.6i'%i)
    stars=read_set_from_file('stars-%6.6i'%i)

    r=(gas.x**2+gas.y**2+gas.z**2)**0.5
    ra,dens=radial_density(r,gas.mass,500)
    ra=numpy.array(map(lambda x: x.number,ra))
    dens=numpy.array(map(lambda x: x.number,dens))
    rad_dens_gas.append([ra,dens])

    r=(stars.x**2+stars.y**2+stars.z**2)**0.5
    ra,dens=radial_density(r,stars.mass,10)
    ra=numpy.array(map(lambda x: x.number,ra))
    dens=numpy.array(map(lambda x: x.number,dens))
    rad_dens_stars.append([ra,dens])

radial_plot(rad_dens_gas,label='gas')
radial_plot(rad_dens_stars,label='stars')

  
