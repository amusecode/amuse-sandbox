import numpy
from matplotlib import pyplot 

from amuse.units import nbody_system
from amuse.units import units,constants

from amuse.datamodel import Particles

from amuse.community.fi.interface import FiViewer

from amuse.io import read_set_from_file

import logging
#logging.basicConfig(level=logging.DEBUG)


i=-1
inp="n"

Mcloud=10000. | units.MSun
Rcloud=1. | units.parsec
conv = nbody_system.nbody_to_si(Mcloud,Rcloud)

viz=FiViewer(conv,redirection="none")

viz.initialize_code()

viz.start_viewer()


while inp!="q":
  if inp=="n":
    i+=1
  if inp=="p":
    i-=1
  if inp.isdigit():
    i=int(inp)    

  print "loading", i

  gas=read_set_from_file('snapshots/gas-%6.6i'%i,'amuse')
  try:
    sink=read_set_from_file('snapshots/sink-%6.6i'%i,'amuse')
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
  
  print "? (n=next snap, p=previous, #=load #, q=quit)"
  inp=raw_input()

