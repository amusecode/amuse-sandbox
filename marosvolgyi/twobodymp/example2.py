from amuse.legacy import *

from interface import twobodympInterface
from interface import twobodymp
import numpy
import time

if __name__ == '__main__':
    instance = twobodympInterface(redirection='none')
    instance.initialization(50)
    print "Precision set to {0}".format( instance.get_precision().precision)
    result = instance.new_particle([0.5, 0.5], [5.0,5.0],
                                   [0.0, 1.0], [0.0, 0.1], [0.0, -0.1],
                                   [0.0, -0.1], [0.0, 0.1], [0.0,-0.2])
    print "Added particle number {0}".format(result.index_of_the_particle)
    instance.commit_particles()
    instance.viewer(1)

    for i in numpy.arange(0.01, 50.0, 0.03):
        instance.plot(0,0,0)
        instance.evolve_system(i)
        instance.plot(255,255,255)
        vx0 = instance.get_velocity(0).vx_
        vy0 = instance.get_velocity(0).vy_
        vz0 = instance.get_velocity(0).vz_
        vx1 = instance.get_velocity(1).vx_
        vy1 = instance.get_velocity(1).vy_
        vz1 = instance.get_velocity(1).vz_
        vy0 -= 0.00
        vy1 += 0.00
        instance.set_velocity(0, vx0, vy0, vz0)
        instance.set_velocity(1, vx1, vy1, vz1)

    s=raw_input()

    instance.stop()
    
