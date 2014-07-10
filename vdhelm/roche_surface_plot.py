"""
    Plot the actual shape of the Roche lobe in a circular orbit.
"""
import numpy

from amuse.units import units, constants, quantities
from amuse.datamodel import Particles

from matplotlib import pyplot
from amuse.plot import plot

length_unit = units.AU
potential_unit = (units.AU/units.yr)**2

def potential_from_mass(x, y, m, a):
    return constants.G * m / ((x + a)**2 + y **2).sqrt()

def potential_at_point(x, y, m1, m2, a1, a2):
    """
        calculate the potential energy at a point close to a binary.
        x and y are location relative to center of mass
        m1 and m2 are masses
        a1 and a2 are distances to center of mass
    """
    potential_1 = potential_from_mass(x, y, m1, a1)
    potential_2 = potential_from_mass(x, y, m2, -a2)
    orbital_potential = 0.5 * (x**2 + y**2) * constants.G * (m1 + m2) / (a1 + a2)**3

    return potential_1 + potential_2 + orbital_potential

def lagrange_points(m1, m2, a1, a2):
    """ calculate the distance from the center of mass, of the L1 and L2 points of mass 2 """
    y= m2*1.0/m1 # y is the mass ratio

    # z is the ratio between the l1 distance and the central object distance
    # first guess using 3z**3 = y
    z = (y/3.)**(1./3.)

    # iteratively get better using 3z**3 = y(1+z)
    for i in range(100):
        z = (y * (1 + z )/3.)**(1./3.)

    l1_r = (a1 + a2) * z
    L1 = a2 - l1_r
    L2 = a2 + l1_r
    return quantities.as_vector_quantity([L1, L2])

def lagrange_point_values(m1, m2, a1, a2):
    x = lagrange_points(m1, m2, a1, a2)
    y = 0 | length_unit
    lagrange_values = potential_at_point(x, y, m1, m2, a1, a2)
    return lagrange_values

def create_grid(plot_range, center, resolution):
    xmin, ymin = center - plot_range/2.
    xmax, ymax = center + plot_range/2.
    x_range = quantities.linspace(xmin, xmax, resolution)
    y_range = quantities.linspace(ymin, ymax, resolution)
    return quantities.as_vector_quantity(numpy.meshgrid(x_range, y_range))

def orbital_angle(binary):
    x1, x2 = binary.x[:2]
    y1, y2 = binary.y[:2]

    x_diff = (x2-x1).value_in(length_unit)
    y_diff = (y2-y1).value_in(length_unit)
    return numpy.arctan2(y_diff, x_diff)

def corotating(x, y, binary):
    X = x - binary.center_of_mass()[0]
    Y = y - binary.center_of_mass()[1]
    # Use the negative angle, because we want to rotate back.
    angle = -orbital_angle(binary)

    co_X = (X * numpy.cos(angle)) - (Y * numpy.sin(angle))
    co_Y = (X * numpy.sin(angle)) + (Y * numpy.cos(angle))
    return co_X, co_Y

def plot_roche_lobe(binary, plot_range, center=[0,0]|length_unit, resolution=200, color='w'):
    """ plot the Roche lobe of the second particle in the binary """
    m1, m2 = binary.mass
    a1, a2 = (binary.center_of_mass() - binary.position).lengths()

    x, y = create_grid(plot_range, center, resolution)
    X, Y = corotating(x, y, binary)

    potentials = potential_at_point(X, Y, m1, m2, a1, a2)

    levels = lagrange_point_values(m1, m2, a1, a2)

    x = x.value_in(length_unit)
    y = y.value_in(length_unit)
    potentials = potentials.value_in(potential_unit)
    levels = levels.value_in(potential_unit)
    return pyplot.contour(x, y, potentials, levels, colors=color, linestyles=('solid', 'dotted'))
