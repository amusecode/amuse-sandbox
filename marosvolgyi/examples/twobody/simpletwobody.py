from amuse.community import *

import sys
sys.path.append('../../twobodymp')
from interface import twobodympInterface
from interface import twobodymp
import numpy

if __name__ == '__main__':
    instance = twobodympInterface()
    instance.initialization(20)

    result = instance.new_particle([0.99, 0.01], [0.1, 0.1],
                                   [0.0, .2], [0.0, 0.0], [0.0,0.0],
                                   [0.0, 0.0], [0.0, 1.0], [0.0,0.0])
    instance.commit_particles()

    for i in numpy.arange(0.01, 5.0, 0.1):
        instance.evolve_system(i)
        pos = instance.get_position(1)
        vel = instance.get_velocity(1)
        print "{0} {1} {2} {3} {4} {5}".format(
            pos.x_, pos.y_, pos.z_,
            vel.vx_, vel.vy_, vel.vz_)

    instance.stop()
    
