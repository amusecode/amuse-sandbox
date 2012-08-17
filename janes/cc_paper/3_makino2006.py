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

"""
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
"""
labels = ["SHARED", "EXTRAPOLATE", "SF", "CC"]#, "CC_KEPLER"] # Can't use CC_KEPLER with softened gravity...
inttypes = [
    huayno2.Huayno.inttypes.SHARED2,
    huayno2.Huayno.inttypes.EXTRAPOLATE,
    huayno2.Huayno.inttypes.HOLD_DKD,
    huayno2.Huayno.inttypes.CC
]
markers = ['', '', '', '']
colors = ['k', 'k', '0.5', 'r']
linestyles = ['-', '--', '-', '--']
eta = 0.1
indices=range(5)

def defNBodyComputation(index=0):
    np.random.seed(index)
    T = 50.0 | nbu.time
    n_samples = 64.0
    return nbe.NBodyComputation(
        initialConditions = plummer.plummer_initial_conditions(n=100), analyticSolution = None,
        dt = T / n_samples, tFinal = T, ndim = 3, outfName = "plummer_cons%02d" % (index,), 
        storeParticles = True
    )

def compute():
    rw = tp.ResultsWriter() # initialize logging
    for index in indices:
        r = defNBodyComputation(index=index) # initialize the problem    
        for (label, inttype) in zip(labels, inttypes):
            hc = huayno2.Huayno()#redirection='none')
            hc.initialize_code()
            hc.parameters.inttype_parameter = inttype
            hc.parameters.timestep_parameter = eta
            hc.parameters.epsilon_squared = 10E-4 | (nbu.length**2) # same as Makino (2006)
            r.runProblemOnIntegrator(hc, label=label, printProgress = True)#, storeHuaynoStats = True)
            r.saveResultsToDisk(rw) # save all results to disk
    del rw # stop logging

def visualize():
    rr = tp.ResultsReader()
    fig_all = plt.figure()
    fig_all.set_size_inches( (12, 10) )
    ax1 = fig_all.add_subplot(221)
    ax2 = fig_all.add_subplot(222)
    ax3 = fig_all.add_subplot(223)
    ax4 = fig_all.add_subplot(224)
    for index in indices:
        r = defNBodyComputation(index)
        r.loadResultsFromDisk(rr)
        r.plotAxTotalEnergyRelativeError(ax1, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
        ax1.set_xlabel('time /$N$-body units/')    
        r.plotAxMomentum(ax2, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
        ax2.set_xlabel('time /$N$-body units/')    
        r.plotAxCOMPos(ax3, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
        ax3.set_xlabel('time /$N$-body units/')    
        r.plotAxAngularMomentum(ax4, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
        ax4.set_xlabel('time /$N$-body units/')    

    #ax.legend(labels_all, bbox_to_anchor=(2.1, 1.05))
    rr.cdto_expdir(subDir="_fig")
    fig_all.savefig("3_makino2006.pdf");
    rr.cdto_prevdir()

if __name__ == '__main__':
    compute()
    visualize()
    print "Fin."
