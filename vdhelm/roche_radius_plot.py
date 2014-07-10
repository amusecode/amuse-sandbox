"""
    Plot the Roche lobe reproducing various plots from Sepinsky.
"""

import numpy
import math
from pylab import arange
from matplotlib import pyplot, cm
from amuse.plot import *
from amuse.ext.roche_radius import Roche_Orbit, sepinsky_formula

def log_range(minimum, maximum, number):
    log_minimum = math.log(minimum)
    log_maximum = math.log(maximum)
    step = (log_maximum - log_minimum)/number
    return [math.exp(x) for x in numpy.arange(log_minimum, log_maximum, step)]

def A_over_nu(eccs, pdf_name=None):
    # Create the data
    plots = []
    labels = []
    nu_range = arange(0.0, 2.0 * numpy.pi, 0.1)
    orbit = Roche_Orbit()
    for ecc in eccs:
        orbit.eccentricity = ecc
        A_values = []
        for nu in nu_range:
            orbit.true_anomaly = nu
            A_values.append(orbit.A)

        plots += [nu_range, A_values]
        labels.append("e=" + str(ecc))


    # Plot the data
    fig = pyplot.figure(figsize = (5, 5))
    semilogy(*plots)
    xlabel("nu")
    ylabel("A(1, e, nu)")
    pyplot.xlim(0, 2 * numpy.pi)
    pyplot.ylim(0.5, 2e4)
    pyplot.legend(labels)
    if pdf_name:
        pyplot.savefig(pdf_name)
    pyplot.show()


def RL_in_Egg_over_A(q_values, xlim=None, ylim=None, pdf_name=None):
    plots = []
    labels = []
    A = numpy.array(log_range(1e-2, 1e4, 100))
    for q in q_values:

        RL_values =sepinsky_formula(q, A)

        plots += [A, RL_values]
        labels.append("q=" + str(q))

    # Plot the data
    plot_data(plots, labels, "A", "RL/RL_egg", xlim, ylim, pdf_name)

def RL_in_Egg_over_q(A_values, xlim=None, ylim=None, pdf_name=None):
    plots = []
    labels = []
    q = log_range(1e-8, 1e8, 100)

    for A in A_values:
        RL_values = sepinsky_formula(q, A)

        plots += [q, RL_values]
        labels.append("A=" + str(A))

    plot_data(plots, labels, "Mass ratio", "RL/RL_egg", xlim, ylim, pdf_name)

def RL_in_D_over_q(A_values=None, f_values=[1.0], e_values=[0.0], xlim=None, ylim=None, pdf_name=None):
    plots = []
    labels = []
    orbit = Roche_Orbit()
    orbit.mass_1 = log_range(1e-8, 1e8, 100) | units.MSun

    if A_values:
        f_values = numpy.sqrt(A_values)

    for f in f_values:
        for e in e_values:
            orbit.angular_velocity_ratio = f
            orbit.eccentricity = e
            RL_values = orbit.sepinsky_roche_radius() / orbit.separation()

            plots += [orbit.mass_ratio, RL_values]
            labels.append("A=" + str(orbit.A))

    plot_data(plots, labels, "Mass ratio", "RL/D", xlim, ylim, pdf_name)

def plot_data(plots, labels, xlab, ylab, xlim, ylim, pdf_name):
    fig = pyplot.figure(figsize = (10, 10))
    semilogx(*plots)
    xlabel(xlab)
    ylabel(ylab)
    if xlim:
        pyplot.xlim(xlim)
    if ylim:
        pyplot.ylim(ylim)
    pyplot.legend(labels, loc=3)
    if pdf_name:
        pyplot.savefig(pdf_name)
    pyplot.show()


def RL_over_q_and_A(pdf_name=None):
    N = 100
    q_range = log_range(1e-8, 1e8, N)
    A_range = log_range(1e-2, 1e4, N)
    X, Y = numpy.meshgrid(q_range, A_range)

    RL_values = sepinsky_formula(X, Y)


    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')

    X, Y = numpy.meshgrid(numpy.log10(q_range), numpy.log10(A_range))
    ax.plot_surface(X, Y, RL_values, cmap=cm.winter)
    pyplot.xlabel("Mass ratio")
    pyplot.ylabel("A")
    pyplot.show()


if __name__ in ('__main__', '__plot__'):

    # recreate Fig. 3 of Sepinski
    A_over_nu([0.0, 0.1, 0.3, 0.6, 0.9])

    # recreate Fig. 9 of Sepinski
    RL_in_Egg_over_q([0.25, 0.5, 1.0, 2.0, 5.0], xlim=[1e-9, 1e9], ylim=[0.5, 1.2])

    # recreate Fig. 8 of Sepinski
    RL_in_D_over_q(A_values=[0.5, 1.0, 2.0, 5.0, 25.0, 100.0, 1000.0], xlim=[1e-9, 1e9], ylim=[-0.05, 1.02])
    # same plot but setting A by changing eccentricity (numbers are not perfect)
    RL_in_D_over_q(e_values=[0.0, 0.0, 0.0994, 0.22963, 0.44256, 0.59762,  0.78369], xlim=[1e-9, 1e9], ylim=[-0.05, 1.02])

    # RL_over_q_and_A()
