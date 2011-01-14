from amuse.community import *

from interface import twobodympInterface
from interface import twobodymp
import numpy
import time
from mpmath import *

PRECISION = 70

def tidalfield(x,y,z, phi, teta, psi):
    F = matrix([[-2*x], [y], [z]])
    Rotx = matrix([[1, 0, 0],[0,cos(phi), -sin(phi)],[0,sin(phi),cos(phi)]])
    Roty = matrix([[cos(teta), 0, sin(teta)],[0,1, 0],[-sin(teta),0,cos(teta)]])
    Rotz = matrix([[cos(psi), -sin(psi), 0],[sin(psi),cos(psi), 0],[0,0,1]])
    Rot = Rotx*Roty*Rotz
    return Rot*F

if __name__ == '__main__':

    instance = twobodympInterface(redirection='none')
    instance.initialization(PRECISION)
    mp.dps = PRECISION
    result = instance.new_particle([0.5, 0.5], [1,1],
                                   [0.0, 1.0], [0.0, 0], [0.0, 0.0],
                                   [0.0, 0.0], [-0.2, 0.2], [0.0,0.0])
    instance.commit_particles()
    instance.viewer(1)

    for i in numpy.arange(0.01, 50.0, 0.01):
        instance.plot(150,50,50)
        instance.evolve_system(i)
        instance.plot(255,255,255)

        x0mp = "".zfill(120)
        y0mp = "".zfill(120)
        z0mp = "".zfill(120)
        vx0mp = "".zfill(120)
        vy0mp = "".zfill(120)
        vz0mp = "".zfill(120)
        x1mp = "".zfill(120)
        y1mp = "".zfill(120)
        z1mp = "".zfill(120)
        vx1mp = "".zfill(120)
        vy1mp = "".zfill(120)
        vz1mp = "".zfill(120)

        x0mp, y0mp, z0mp, res = instance.get_position_mp(0, x0mp,y0mp,z0mp, PRECISION)
        vx0mp, vy0mp, vz0mp, res = instance.get_velocity_mp(0, vx0mp,vy0mp,vz0mp, PRECISION)
        x1mp, y1mp, z1mp, res = instance.get_position_mp(1, x1mp,y1mp,z1mp, PRECISION)
        vx1mp, vy1mp, vz1mp, res = instance.get_velocity_mp(1, vx1mp,vy1mp,vz1mp, PRECISION)

        vx0 = mpf(vx0mp[0])#instance.get_velocity(0).vx_
        vy0 = mpf(vy0mp[0])#instance.get_velocity(0).vy_
        vz0 = mpf(vz0mp[0])#instance.get_velocity(0).vz_
        vx1 = mpf(vx1mp[0])#instance.get_velocity(0).vx_
        vy1 = mpf(vy1mp[0])#instance.get_velocity(0).vy_
        vz1 = mpf(vz1mp[0])#instance.get_velocity(0).vz_

        x0 = mpf(x0mp[0])
        y0 = mpf(y0mp[0])
        z0 = mpf(z0mp[0])
        x1 = mpf(x1mp[0])
        y1 = mpf(y1mp[0])
        z1 = mpf(z1mp[0])

        scale = 0.001

        F1 = scale * tidalfield(x0,y0,z0, 0, 0, 0)
        vx0 += F1[0]
        vy0 += F1[1]
        vz0 += F1[2]
        F1 = scale * tidalfield(x1,y1,z1, 0, 0, 0)
        vx1 += F1[0]
        vy1 += F1[1]
        vz1 += F1[2]

        instance.set_velocity_mp(0, nstr(vx0,100), nstr(vy0,100), nstr(vz0,100))
        instance.set_velocity_mp(1, nstr(vx1,100), nstr(vy1,100), nstr(vz1,100))
        
        print x0,y0,z0,x1,y1,z1
        
    s=raw_input()

    instance.stop()
    
