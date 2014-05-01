import Image
import numpy

from matplotlib import pyplot

from amuse.datamodel import Particles
from amuse.units import units,constants

from amuse.ic.salpeter import new_salpeter_mass_distribution

from ubvivariant import rgb_frame

from amuse.community.sse.interface import SSE

def generate_set(Ntarget=1000):
  a=Image.open("amuse.png")
  
  Nx,Ny=a.size
  
  a=numpy.array(a.getdata())
  
  print a.max()
  
  b=(a[:,0]**2/3 + a[:,1]**2/3 + a[:,2]**2/3)**0.5
  
  b=b.reshape((Ny,Nx))
  
  
  numpy.random.seed(123332)
  
  parts=Particles()
  
  N=0
  while len(parts)<Ntarget:
    print len(parts)
    Nsample=3*(Ntarget-N)/2
    x=numpy.random.uniform(0.,Nx,Nsample)
    y=numpy.random.uniform(0.,Ny,Nsample)
    ix=numpy.int64(x)
    iy=numpy.int64(y)
    a=numpy.where( b[iy,ix]<128)[0][0:Ntarget-len(parts)]
    new=Particles(len(a),x=x[a],y=-y[a])
    parts.add_particles(new)
    
    
  stars=Particles(len(parts))
  stars.x=-(parts.x-Nx/2.) | units.parsec
  stars.z=-(parts.y+Ny/2.) | units.parsec
  stars.mass=new_salpeter_mass_distribution(
        len(parts),
        mass_min = 0.3 | units.MSun,
        mass_max = 80.0 | units.MSun,
        alpha = -2.35#,random=notsorandom()
    )
  stars.y=0.*stars.x
  return stars

def evolve(parts,tend=0| units.Myr):
  sse=SSE()
  parts_in_code=sse.particles.add_particles(parts)
  xmin=parts.x.min()
  xmax=parts.x.max()
  for i,p in enumerate(parts_in_code):
    p.evolve_for( tend* (xmax-parts[i].x)/(xmax-xmin) )
    print i,p.age

#  sse.evolve_model(tend)
  channel=sse.particles.new_channel_to(parts)
  channel.copy_attributes(["luminosity","radius"])

def calculate_effective_temperature(luminosity, radius):
    return ((luminosity/(constants.four_pi_stefan_boltzmann*radius**2))**.25).in_(units.K)

if __name__=="__main__":
  stars=generate_set(20000)
  evolve(stars, tend=100.| units.Myr)
  
  stars.temperature=calculate_effective_temperature(stars.luminosity,stars.radius)  
  
  label="amuse"
  L=800. | units.parsec
  image_size=[960*2,2*540]
  percentile=0.998
  
#  vmax=rgb_frame(stars,dryrun=True,image_width=L,image_size=image_size,percentile=percentile)

  vmax,image=rgb_frame(stars,dryrun=False, image_width=L, image_size=image_size, 
    percentile=percentile,multi_psf=False)  
  
  imrgb=Image.fromstring(image['mode'], image['size'], image['pixels'])
  imrgb.save(label+"-rgb.png")
    
  
  
    
