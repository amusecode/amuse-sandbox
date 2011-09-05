"""
Amuse-implementation of complex N-body problems with a known analytic solution as described in:
    [MN08] Cristopher Moore and Michael Nauenberg. New Periodic Orbits for the n-Body Problem. 
    URL: http://arxiv.org/abs/math/0511219v2

Implemented problems:
1) Cuboctahedron orbit only for the special case of m=1.
2) Planar criss-cross orbit from.

"""

import math
import numpy as np

from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE
from .. nbody_experiments import *


def cuboctahedron_initial_conditions(m = 1):
    """
    Implementation of the cuboctahedron-shaped orbits from Moore & Nauenberg, 2008, Section 2
    """

    x = [[0, -0.69548, 0.69548]]
    x.append([+x[0][0], -x[0][1], -x[0][2]])
    x.append([-x[0][0], +x[0][1], -x[0][2]])
    x.append([-x[0][0], -x[0][1], +x[0][2]])

    v = [[0.87546, -0.31950, -0.31950]]
    v.append([+v[0][0], -v[0][1], -v[0][2]])
    v.append([-v[0][0], +v[0][1], -v[0][2]])
    v.append([-v[0][0], -v[0][1], +v[0][2]])

    for t_shift in map(lambda j: 2 * j * math.pi / m, range(1, m)):

        t_shift = t_shift
        x0 = cuboctahedron_f(t_shift)
        y0 = cuboctahedron_f(t_shift + 2 * math.pi / 3)
        z0 = cuboctahedron_f(t_shift + 4 * math.pi / 3)

        vx = cuboctahedron_f(t_shift)
        vy = cuboctahedron_f(t_shift + 2 * math.pi / 3)
        vz = cuboctahedron_f(t_shift + 4 * math.pi / 3)
    
        x.append([+x0, +y0, +z0])
        x.append([+x0, -y0, -z0])
        x.append([-x0, +y0, -z0])
        x.append([-x0, -y0, +z0])
        
        v.append([+vx, +vy, +vz])
        v.append([+vx, -vy, -vz])
        v.append([-vx, +vy, -vz])
        v.append([-vx, -vy, +vz])

    stars = data.Stars(len(x))
    for i in range(len(stars)):
        star = stars[i]
        star.mass = 1. | nbody_system.mass
        star.position = x[i] | nbody_system.length
        star.velocity = v[i] | ( nbody_system.length / nbody_system.time )
        star.radius = 0. | nbody_system.length # obligatory for Hermite() integrator

    return stars

"""
Fourier' coefficients for the analytical solutions of the orbits
key determines the number of particles on a single orbit
"""
cuboctahedron_a_all = {
     1:   [1.00000, 0.03282, -.00098, -.00036, -.00003],
     3:   [1.00000, -.04629, -.00472, -.00269, 0.00056, 0.00010, 0.00007, -.00002],
     5:   [1.00000, -.03991, 0.00359, 0.00161, -.00168, -.00043, -.00028, 0.00012, 0.00001, 0.00001],
     7:   [1.00000, -.03335, 0.00885, 0.00445, -.00334, -.00039, -.00028, -.00013, -.00011, -.00008, -.00005],
     61:  [1.00000, -.02302, 0.01054, 0.00544, -.00181],
     81:  [1.00000, -.02285, 0.01025, 0.00522, -.00175],
     101: [1.00000, -.02265, 0.01008, 0.00507, -.00167],
     121: [1.00000, -.02235, 0.00991, 0.00494, -.00156],
}

def cuboctahedron_f(t, m = 1, normalize = True):

    a = cuboctahedron_a_all[m]
    
    f = 0
    for i in range((len(a))):
        k = 2 * i + 1
        f += a[i] * math.sin(k * t)

    # magic for un-normalizing the Fourier' coefficients given in the article
    if normalize: 
        C = -0.69548 / cuboctahedron_f(2.0 / 3.0 * math.pi, normalize = False)
        f = C * f

    return f

def cuboctahedron_fdot(t, m = 1, normalize = True):

    a = cuboctahedron_a_all[m]
    
    f = 0
    for i in range((len(a))):
        k = 2 * i + 1
        f += a[i] * k * math.cos(k * t)
        
    # magic for un-normalizing the Fourier' coefficients given in the article
    if normalize: 
        C = -0.69548 / cuboctahedron_f(2.0 / 3.0 * math.pi, normalize = False)
        f = C * f

    return f

def cuboctahedron_analytic_orbits(t, m = 1):

    # flipping the sign of t magically makes the analytic trajectory match the integrated trajectory
    # in time; no clue as to why it works
    t = -np.array(t)

    x = np.array(map(cuboctahedron_f, t))
    y = np.array(map(cuboctahedron_f, t + 2 * math.pi / 3))
    z = np.array(map(cuboctahedron_f, t + 4 * math.pi / 3))
    
    return [
            [ x,  y,  z],
            [+x, -y, -z],
            [-x, +y, -z],
            [-x, -y, +z],
            ]

def planar_crisscross_initial_conditions():
    """
    Implementation of the planar criss-cross orbit from Moore & Nauenberg, 2008, Section 3
    There is numerical proof that the orbit is stable.
    """

    x = [ 
         [+1.07590, 0, 0], 
         [-0.07095, 0, 0], 
         [-1.00496, 0, 0], 
    ]

    v = [ 
         [0, +0.19509, 0], 
         [0, -1.23187, 0], 
         [0, +1.03678, 0], 
    ]
    
    m = [1, 1, 1]
    
    return particles_from_floats(m, x, v)

    stars = data.Stars(len(x))
    for i in range(len(stars)):
        star = stars[i]
        star.mass = 1.0 | nbody_system.mass
        star.position = x[i] | nbody_system.length
        star.velocity = v[i] | ( nbody_system.length / nbody_system.time )
        star.radius = 0.0 | nbody_system.length # obligatory for Hermite() integrator

    return stars
 
def planar_crisscross_analytic_solution(t):

    a1 = [1.09764, -.02809, 0.00724, -.00121, 0.00040, -.00010, 0.00003, -.00001]
    a3 = [-.98868, -.00442, -.01100, -.00010, -.00069, -.00001, -.00006, 0.00000]
    b1 = [0.10896, 0.03251, -.00376, 0.00131, -.00029, 0.00010, -.00003, 0.00001]

    a2 = [-b1[0], +b1[1], -b1[2], +b1[3], -b1[4], +b1[5], -b1[6], -b1[7]]
    b2 = [-a1[0], +a1[1], -a1[2], +a1[3], -a1[4], +a1[5], -a1[6], -a1[7]]
    b3 = [-a3[0], +a3[1], -a3[2], +a3[3], -a3[4], +a3[5], -a3[6], -a3[7]]

    def _fcos(t, c):
        """
        essentially, implementation of formula (11) from the paper
        """
        
        f = 0
        for i in range((len(c))):
            k = 2 * i + 1
            f += c[i] * math.cos(k * t)
            
        return f

    def _fsin(t, c):
        """
        essentially, implementation of formula (11) from the paper
        """
        
        f = 0
        for i in range((len(c))):
            k = 2 * i + 1
            f += c[i] * math.sin(k * t)
            
        return f

    x1 = np.array(map(lambda t: _fcos(t, a1), t))
    y1 = np.array(map(lambda t: _fsin(t, b1), t))
    x2 = np.array(map(lambda t: _fcos(t, a2), t))
    y2 = np.array(map(lambda t: _fsin(t, b2), t))
    x3 = np.array(map(lambda t: _fcos(t, a3), t))
    y3 = np.array(map(lambda t: _fsin(t, b3), t))
    z = np.zeros(len(x1))

    return [[x1, y1, z], [x2, y2, z], [x3, y3, z]]


if __name__ == '__main__':

    # cuboctahedron, m=1, 1/3 orbit
    r = NBodyComputation(
        initialConditions = cuboctahedron_initial_conditions(), analyticSolution = cuboctahedron_analytic_orbits,
        dt = 0.1 | nbody_system.time, tFinal = 2 | nbody_system.time,
        ndim = 3, outfName = "MN08_cuboctahedron_T0.5"
    )
    r.runProblemOnIntegrator(PhiGRAPE(), label="PhiGRAPE", printProgress = False)
    r.runProblemOnIntegrator(Hermite(), label="Hermite", printProgress = False)

    r.plotTotalEnergyRelativeError(labels = ["PhiGRAPE", "Hermite"])
    r.plotTrajectories(label="PhiGRAPE")
    r.plotTrajectories(label="Hermite")

    # cuboctahedron, m=1, full orbit
    r = NBodyComputation(
        initialConditions = cuboctahedron_initial_conditions(), analyticSolution = cuboctahedron_analytic_orbits,
        dt = 0.1 | nbody_system.time, tFinal = 6.2 | nbody_system.time,
        ndim = 3, outfName = "MN08_cuboctahedron_T1.0"
    )
    r.runProblemOnIntegrator(PhiGRAPE(), label="PhiGRAPE", printProgress = False)
    r.runProblemOnIntegrator(Hermite(), label="Hermite", printProgress = False)

    r.plotTotalEnergyRelativeError(labels = ["PhiGRAPE", "Hermite"])
    r.plotTrajectories(label="PhiGRAPE")
    r.plotTrajectories(label="Hermite")


    # planar criss-cross, 1/2 orbit
    r = NBodyComputation(
        initialConditions = planar_crisscross_initial_conditions(), analyticSolution = planar_crisscross_analytic_solution,
        dt = 0.1 | nbody_system.time, tFinal = 3.1 | nbody_system.time,
        ndim = 2, outfName = "MN08_crisscross_T0.5"
    )
    r.runProblemOnIntegrator(PhiGRAPE(), label="PhiGRAPE", printProgress = False)
    r.runProblemOnIntegrator(Hermite(), label="Hermite", printProgress = False)

    r.plotTotalEnergyRelativeError(labels = ["PhiGRAPE", "Hermite"])
    r.plotTrajectories(label="PhiGRAPE")
    r.plotTrajectories(label="Hermite")

    # planar criss-cross, full orbit
    r = NBodyComputation(
        initialConditions = planar_crisscross_initial_conditions(), analyticSolution = planar_crisscross_analytic_solution,
        dt = 0.1 | nbody_system.time, tFinal = 6.2 | nbody_system.time,
        ndim = 2, outfName = "MN08_crisscross_T1.0"
    )
    r.runProblemOnIntegrator(PhiGRAPE(), label="PhiGRAPE", printProgress = False)
    r.runProblemOnIntegrator(Hermite(), label="Hermite", printProgress = False)

    r.plotTotalEnergyRelativeError(labels = ["PhiGRAPE", "Hermite"])
    r.plotTrajectories(label="PhiGRAPE")
    r.plotTrajectories(label="Hermite")
