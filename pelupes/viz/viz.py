import numpy
from matplotlib import pyplot 

from amuse.units import nbody_system
from amuse.units import units,constants

from amuse.datamodel import Particles

from amuse.community.fi.interface import FiViewer,FiMap

from amuse.io import read_set_from_file

import logging
#logging.basicConfig(level=logging.DEBUG)

def column_density_map(convert,part, image_target,viewpoint,image_angle, ratio, N=500, projection="perspective"):  
  mapper=FiMap(convert)
  
  mapper.parameters.minimum_distance=1. | units.AU
  mapper.parameters.image_size=[N,int(N/ratio)]
  mapper.parameters.image_target=image_target

  image_width=2*(image_target-viewpoint).length()*numpy.tan(numpy.pi*image_angle/180./2)
  mapper.parameters.image_width=image_width
  mapper.parameters.projection_direction=(image_target-viewpoint)/(image_target-viewpoint).length()
  mapper.parameters.projection_mode=projection
  mapper.parameters.image_angle=image_angle
  mapper.parameters.viewpoint=viewpoint

  print mapper.parameters

  mapper.commit_parameters()

  part.weight=part.mass.value_in(units.amu)

  mapper.particles.add_particles(part)
    
  projected=mapper.image.pixel_value 
        
  mapper.stop()

  im=numpy.transpose(projected)
  
  return im, image_width

def makemap(convert,viz,projection="perspective"):
  image_target=viz.parameters.image_target
  viewpoint=viz.parameters.viewpoint
  angle=viz.parameters.image_angle
  ratio=viz.parameters.image_ratio
  
  print "ratio, angle:", ratio,angle
  
  xangle=2*180/numpy.pi*numpy.arctan(ratio*numpy.tan(numpy.pi*angle/180/2))
  print "xangle:",xangle
  
  im,image_width=column_density_map(convert,viz.gas_particles.copy(),image_target,viewpoint,xangle,ratio,
        projection=projection)

  iw=image_width.value_in(units.AU)

  extent=[-iw,iw,-iw/ratio,iw/ratio]

  md=numpy.mean(im)
  im=numpy.log10(im+md/1.e6)
  
  vmax=numpy.log10(md)+2
  vmin=numpy.log10(md)-2
  
  pyplot.clf()
  pyplot.imshow(im, vmin=vmin,vmax=vmax,origin="lower",cmap="copper",extent=extent)
  pyplot.draw()
  

def update_view(viz,i):    
  print "loading", i

  gas=read_set_from_file(snapdir+'/gas-%6.6i'%i,'amuse')
  try:
    sink=read_set_from_file(snapdir+'/sink-%6.6i'%i,'amuse')
  except:
    sink=Particles()

  if len(viz.gas_particles)>0:
    viz.gas_particles.remove_particles(viz.gas_particles)
  if len(viz.dm_particles)>0:
    viz.dm_particles.remove_particles(viz.dm_particles)

  viz.gas_particles.add_particles(gas)
  if len(sink):
    viz.dm_particles.add_particles(sink)

#  viz.trigger_viewer_refresh()


i=0
inp="n"

Mcloud=50. | units.MSun
Rcloud=0.5*0.375 | units.parsec
conv = nbody_system.nbody_to_si(Mcloud,Rcloud)

viz=FiViewer(conv,redirection="none")

viz.initialize_code()

viz.start_viewer()

snapdir="snapshots"

f=pyplot.figure(figsize=(8,6))

pyplot.ion()
pyplot.show()

iprev=-1

while inp!="q":

  do_makemap=False
  if inp[-1:] in ["m","M"]:
     do_makemap=True
     projection="perspective" if inp[-1:]=="m" else "parallel" 
     inp=inp[:-1]

  if inp=="n":
    i+=1
  if inp=="p":
    i-=1
  if inp.isdigit():
    i=int(inp)    

  if i!=iprev:
    iprev=i
    update_view(viz,i)
  if do_makemap:
    makemap(conv,viz,projection=projection)
  
  print "? (n=next snap, p=previous, #=load #,m=map, q=quit)"
  inp=raw_input()
