from amuse.units import units
from amuse.units import nbody_system
import numpy as np

from amuse.support import data
def Pythagorean_initial_conditions(): 
    """ 
    Initial conditions for the Pythagorean system, from:
        amuse/sandbox/spz/Pythagorean.py
    """

    stars = data.Stars(3)
        
    unit_velocity = nbody_system.length/nbody_system.time
    star0 = stars[0]
    star0.mass = 3.0 | nbody_system.mass
    star0.position = [1, 3, 0] | nbody_system.length
    star0.velocity = [0, 0, 0] | unit_velocity
    star0.radius = 0.0 | nbody_system.length

    star1 = stars[1]
    star1.mass = 4.0 | nbody_system.mass
    star1.position = [-2, -1, 0] | nbody_system.length
    star1.velocity = [0, 0, 0] | unit_velocity
    star1.radius = 0.0 | nbody_system.length

    star2 = stars[2]
    star2.mass = 5.0 | nbody_system.mass
    star2.position = [1, -1, 0] | nbody_system.length
    star2.velocity = [0, 0, 0] | unit_velocity
    star2.radius = 0.0 | nbody_system.length

    return stars
