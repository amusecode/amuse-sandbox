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

labels =  ["EXTRAPOLATE", "SF", "CC"]
inttypes = [
    huayno2.Huayno.inttypes.EXTRAPOLATE,
    huayno2.Huayno.inttypes.HOLD_DKD,
    huayno2.Huayno.inttypes.CC
]
markers = ['', '', '']
colors = ['k', '0.5', 'r']
linestyles = ['-', '-', '--']

eta=0.1; no_particles=20; T=700.0|nbu.time; n_samples=512.0 # test run with 100 bodies
#eta=0.01; no_particles=1024; T=700.0|nbu.time; n_samples=512.0

def defNBodyComputation(index=0):
    np.random.seed(index)
    return nbe.NBodyComputation(
        initialConditions = plummer.plummer_initial_conditions(n=no_particles), analyticSolution = None,
        dt = T / n_samples, tFinal = T, ndim = 3, outfName = "plummer_cons%02d" % (index,), 
        storeParticles = True
    )

def compute():
    rw = tp.ResultsWriter() # initialize logging
    r = defNBodyComputation() # initialize the problem    
    for (label, inttype) in zip(labels, inttypes):
        hc = huayno2.Huayno()#redirection='none')
        hc.initialize_code()
        hc.parameters.inttype_parameter = inttype
        hc.parameters.timestep_parameter = eta
        hc.parameters.epsilon_squared = (1/256.0)**2 | (nbu.length**2) # same as Makino (2006)
        r.runProblemOnIntegrator(hc, label=label, printProgress=True, storeHuaynoStats=True)
        r.saveResultsToDisk(rw) # save all results to disk
    del rw # stop logging

def visualize():
    rr = tp.ResultsReader()
    fig_all = plt.figure()
    fig_all.set_size_inches( (12, 18) )
    r = defNBodyComputation()
    r.loadResultsFromDisk(rr)

    print "Plotting energy..."
    ax1 = fig_all.add_subplot(521)
    r.plotAxTotalEnergyRelativeError(ax1, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
    ax1.set_xlabel('time /$N$-body units/')    

    print "Plotting momentum..."
    ax2 = fig_all.add_subplot(522)
    r.plotAxMomentum(ax2, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
    ax2.set_xlabel('time /$N$-body units/')    

    print "Plotting COM..."
    ax3 = fig_all.add_subplot(523)
    r.plotAxCOMPos(ax3, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
    ax3.set_xlabel('time /$N$-body units/')    

    print "Plotting angular momentum..."
    ax4 = fig_all.add_subplot(524)
    r.plotAxAngularMomentum(ax4, labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
    ax4.set_xlabel('time /$N$-body units/')    

    print "Plotting radii..."
    ax5 = fig_all.add_subplot(525)
    #r.plotAxPlummerRadii(ax5, labels=labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
    #r.plotAxCoreRadius(ax5, labels=labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
    ax5.set_xlabel('time /$N$-body units/')    
    ax5.set_ylabel('{Lagrangian, core} radius')

    print "Plotting density..."
    ax6 = fig_all.add_subplot(526)
    #r.plotAxCoreDensity(ax6, labels=labels, logYScale=True, colors=colors, markers=markers, linestyles=linestyles)
    ax6.set_xlabel('time /$N$-body units/')
    ax6.set_ylabel('core density')

    print "Plotting kicks..."
    ax7 = fig_all.add_subplot(527)
    ax7.set_xlabel('time /$N$-body units/')
    ax7.set_ylabel('no of kick evaluations')
    r.plotAxTotalKicks(ax7, labels=labels, logYScale=False, colors=colors, markers=markers, linestyles=linestyles)

    print "Plotting drifts..."
    ax8 = fig_all.add_subplot(528)
    ax8.set_xlabel('time /$N$-body units/')
    ax8.set_ylabel('no of drift evaluations')
    r.plotAxTotalDrifts(ax8, labels=labels, logYScale=False, colors=colors, markers=markers, linestyles=linestyles)

    print "Plotting TS..."
    ax9 = fig_all.add_subplot(529)
    ax9.set_xlabel('time /$N$-body units/')
    ax9.set_ylabel('no of time step evaluations')
    r.plotAxTotalTS(ax9, labels=labels, logYScale=False, colors=colors, markers=markers, linestyles=linestyles)

    print "Plotting wall-clock..."
    ax0 = fig_all.add_subplot(5,2,10)
    ax0.set_xlabel('time /$N$-body units/')
    ax0.set_ylabel('wall-clock time')
    r.plotAxWallClockTime(ax0, labels=labels, logYScale=False, colors=colors, markers=markers, linestyles=linestyles)

    #ax.legend(labels_all, bbox_to_anchor=(2.1, 1.05))
    rr.cdto_expdir(subDir="_fig")
    fig_all.savefig("4_nitadori2008.pdf");
    rr.cdto_prevdir()

if __name__ == '__main__':
    #compute()
    visualize()
    print "Fin."

