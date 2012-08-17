#!/usr/bin/env amuse.sh

import math
import sys
import time
import numpy as np
import amuse.support.units.nbody_system as nbu
from amuse.support.units import units
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.ph4.interface import ph4
import shared_code.nbody_experiments.nbody_experiments as nbe
import shared_code.nbody_experiments.problems.two_body as tbp
import shared_code.nbody_experiments.problems.plummer_with_planets as plummer
import shared_code.nbody_experiments.problems.mn08_problems as mn08
import shared_code.nbody_experiments.problems.pythagorean_test as pt
import amuse.community.huayno.interface as huayno2
import shared_code.ductape as tp
import matplotlib.pyplot as plt

labels = ["SHARED", "EXTRAPOLATE", "SF", "CC", "CC_KEPLER"]
inttypes = [
    huayno2.Huayno.inttypes.SHARED2,
    huayno2.Huayno.inttypes.EXTRAPOLATE,
    huayno2.Huayno.inttypes.HOLD_DKD,
    huayno2.Huayno.inttypes.CC,
    huayno2.Huayno.inttypes.CC_KEPLER
]
markers = ['', '', '', '', '']
colors = ['k', 'k', '0.5', 'r', 'r']
linestyles = ['-', '--', '-', '--', '-.']
eta = 0.01

def defNBodyComputation():
    T = 10*2*math.pi | nbu.time
    n_samples = 84
    return nbe.NBodyComputation(
        initialConditions=mn08.planar_crisscross_initial_conditions(), 
        analyticSolution=mn08.planar_crisscross_analytic_solution,
        dt = T / float(n_samples), tFinal = T, ndim = 2, outfName = "crisscross_cons"
    )

def compute():
    rw = tp.ResultsWriter() # initialize logging
    r = defNBodyComputation() # initialize the problem
    # Huayno integrators
    for (label, inttype) in zip(labels, inttypes):
        hc = huayno2.Huayno()#redirection='none')
        hc.initialize_code()
        hc.parameters.inttype_parameter = inttype
        hc.parameters.timestep_parameter = eta
        r.runProblemOnIntegrator(hc, label=label, printProgress = True)#, storeHuaynoStats = True)
    r.saveResultsToDisk(rw) # save all results to disk
    del rw # stop logging

def visualize():
    xticks = np.arange(0, 11*2*np.pi, 2*2*np.pi)
    def _orbitalPeriodFormatter(x, pos):
        return "%.0fT" % (x / (2*np.pi))
    orbitalPeriodFormatter = plt.FuncFormatter(_orbitalPeriodFormatter)
  
    # load results
    rr = tp.ResultsReader()
    r = defNBodyComputation()
    r.loadResultsFromDisk(rr)

    # plot a comparison in energy
    rr.cdto_expdir(subDir="_fig")

    # conserved quantities
    fig = plt.figure()
    fig.set_size_inches( (12, 12) )
    ax = fig.add_subplot(321)
    ax.set_xlim(0, 10)
    #ax.set_ylim((10E-18, 10E-5))
    r.plotAxTotalEnergyRelativeError(ax, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_major_formatter(orbitalPeriodFormatter)

    ax = fig.add_subplot(322)
    ax.set_xlim(0, 10)
    ax.set_ylim((10E-19, 10E-6))
    r.plotAxMomentum(ax, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_major_formatter(orbitalPeriodFormatter)

    ax = fig.add_subplot(323)
    ax.set_xlim(0, 10)
    #ax.set_ylim((-10E-18, +10E-18))
    r.plotAxCOMVel(ax, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_major_formatter(orbitalPeriodFormatter)

    ax = fig.add_subplot(324)
    ax.set_xlim(0, 10)
    ax.set_ylim((10E-18, 10E-6))
    r.plotAxAngularMomentum(ax, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
    ax.set_xlabel('integration time $t$')    
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_major_formatter(orbitalPeriodFormatter)

    ax = fig.add_subplot(325)
    ax.set_xlim(0, 10)
    r.plotAxTrajError(ax, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
    ax.set_xlabel('integration time $t$')
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_major_formatter(orbitalPeriodFormatter)

    ax.legend(labels, bbox_to_anchor=(1.6, 1.0))
    fig.savefig("2b_crisscross_cons.pdf")
    rr.cdto_prevdir()

if __name__ == '__main__':
    compute()
    visualize()
    print "Fin."
