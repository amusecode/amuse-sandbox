import sys
from pylab import *

from amuse.community.bhtree.interface import BHTreeInterface, BHTree
from amuse.community.hermite0.interface import HermiteInterface, Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPEInterface, PhiGRAPE

from amuse.units import nbody_system
from amuse.units import units

from amuse.support.codes.core import is_mpd_running

from amuse import datamodel
def Pythagorean_system() : 

    stars = datamodel.Stars(3)
        
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
        
if __name__ == '__main__':

    assert is_mpd_running()

    stars = Pythagorean_system()

    end_time = 40 | nbody_system.time
    dt = 0.0625 | nbody_system.time 
    convert_nbody = None

    gravity = Hermite()
#    gravity = BHTree()
    if isinstance(gravity, BHTree) :
        gravity.initialize_code()
        eps = 0.001 | nbody_system.length
    else :
        eps = 0 | nbody_system.length
    gravity.parameters.epsilon_squared = eps**2 

    stars.radius = 0.0 | nbody_system.length

    gravity.particles.add_particles(stars)
    gravity.commit_particles()
    from_model_to_gravity = stars.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(stars)
    
    time = 0.0 | nbody_system.time
    Ek0 = gravity.kinetic_energy
    Ep0 = gravity.potential_energy
    Et0 = Ek0 + Ep0
    print "Time= ", time, " E= ", Et0, Ek0, Ep0, " Q= ", Ek0/Ep0

    print "evolving the model until t = " + str(end_time)
    x = [None]*(len(stars))
    y = [None]*(len(stars)) 
    xp = []
    yp = []
    Tlabel = 'x'
    Rlabel = 'y'
    i=0
    while time < end_time:
        i+=1
        time += dt
        
        gravity.evolve_model(time)
        from_gravity_to_model.copy()

        Ek = gravity.kinetic_energy
        Ep = gravity.potential_energy
        Et = Ek + Ep
        print "Time= ", time, " E= ", Et, Ek, Ep, Ek/Ep, (Et-Et0)/Et0
        for si in range(len(stars)) :
            x[si] = stars[si].x.value_in(nbody_system.length)
            y[si] = stars[si].y.value_in(nbody_system.length)
            xp.append(x[si])
            yp.append(y[si])
            if len(xp) > 30 :
                xp.pop(0)
                yp.pop(0)
        print x, y
        errorbarsize = 0.0
        plt.errorbar(xp, yp, errorbarsize, fmt='.') 
        plt.errorbar(x, y, errorbarsize, fmt='o') 
        plt.axis([-5, 5, -5, 5])
        plt.xlabel(Tlabel)
        plt.ylabel(Rlabel)
        filename = str('%03d' % i) + '.png'
        plt.savefig(filename, dpi=100)
        plt.clf()

# Make a movie using
# %> mencoder "mf://*.png" -mf fps=30 -ovc raw -o output.avi
# run it with
# %> mpegplay output.avi
