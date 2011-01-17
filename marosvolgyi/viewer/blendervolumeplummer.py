try:
    from OpenGL.GL import *
    from OpenGL.GLU import *
    import pygame
    from pygame.locals import *
except ImportError:
    if __name__ == '__main__':
        raise Exception("OpenGL and pygame bindings are not installed, cannot run")
    
        
#from OpenGL.GLUT import *

import struct

from amuse.support.data import core
from amuse.support.units import units
from amuse.support.units import nbody_system 
from amuse.ext.plummer import MakePlummerModel
from amuse.ext.kingmodel import MakeKingModel
from amuse.community.hermite0.interface import Hermite
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.octgrav.interface import Octgrav

from datetime import date, timedelta
import time

import sys
import numpy as np
import random

if __name__ == "__main__":

    nstars = int(sys.argv[1])
    workers = int(sys.argv[2])
    method = sys.argv[3]
    print nstars
    seed = None
    stars = MakePlummerModel(nstars).result#, convert_nbody, random_state = seed).result

    if method == 'octgrav':
        gravity = Octgrav()
    elif method == 'phigrape':					
        gravity = PhiGRAPE(PhiGRAPE.NBODY)
    elif method == 'bhtree':
        gravity = BHTree(number_of_workes = workers) 
    elif method == 'hermite':
        gravity = Hermite(
                          number_of_workers = workers
                          #debugger = "xterm",
                          #redirection = "none"
                          )
    gravity.initialize_code()
    gravity.parameters.epsilon_squared = 0.001 | nbody_system.length ** 2
    
    stars.radius = 0.000001 | nbody_system.length
    gravity.particles.add_particles(stars)
    gravity.stopping_conditions.pair_detection.enable()
    from_model_to_gravity = stars.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(stars)
    
    gravity.commit_particles()
    model_time = 0.0
    pairs = 0
    image_count = 0
    
    model_time +=0.1

    gravity.evolve_model(model_time|nbody_system.time)
    from_gravity_to_model.copy()

    #print gravity.kinetic_energy
    #print gravity.potential_energy
    files = open('field.bin', 'wb')

    for xx in np.arange(-1,1,0.05):
        for yy in np.arange(-1,1,0.05):
            for zz in np.arange(-1,1,0.05):
                phi = gravity.get_potential_at_point(0.01|nbody_system.length,
                                                     xx|nbody_system.length, 
                                                     yy|nbody_system.length, 
                                                     zz|nbody_system.length)
                    #print -phi.value_in(nbody_system.length**2/(nbody_system.time**2))
                value = round(40*(-phi.value_in(nbody_system.length**2/(nbody_system.time**2))))
                if value>255: value = 255
                files.write(struct.pack('B', value))
                
    files.flush()
    files.close()


    gravity.stop()
