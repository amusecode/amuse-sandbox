import os
import sys
import numpy

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    
from amuse.units import nbody_system
from amuse.units import units

    
from amuse.community.fi import interface as interface
from amuse.community.phiGRAPE.interface import PhiGRAPE

from amuse.ic.plummer import new_plummer_sphere

number_of_stars=1000

particles = new_plummer_sphere(number_of_stars)
   
gravity = PhiGRAPE(channel_type="ibis",hostname="paddegat")
gravity.parameters.epsilon_squared = 0.01 | nbody_system.length ** 2
    
gravity.particles.add_particles(particles)

times=nbody_system.time([0.,0.2,0.4,0.6,0.8,1.0,1.2,1.4])
f=pyplot.figure(figsize=(8,8))

for i,t in enumerate(times):
    gravity.evolve_model(t)

    x=gravity.particles.x.value_in(nbody_system.length)
    y=gravity.particles.y.value_in(nbody_system.length)

    pyplot.plot(x,y,'r .')
    pyplot.xlim(-3,3)
    pyplot.ylim(-3,3)
    if i==7:
        pyplot.plot(x,y,'b +')
        pyplot.xlabel('nbody length')
            
#    pyplot.show()
pyplot.savefig('cluster.png')


