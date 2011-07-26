import math
import optparse
import sys
import time
import numpy as np
from amuse.support.units import units
import amuse.support.units.nbody_system as nbu
import nbody_experiments.nbody_experiments as nbe
import nbody_experiments.problems.two_body as tbp
import nbody_experiments.problems.plummer_with_planets as plummer
import nbody_experiments.problems.mn08_problems as mn08
import nbody_experiments.problems.pythagorean_test as pt
import ceutils.ceutils as ce
import matplotlib.pyplot as plt
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.ph4.interface import ph4
from amuse.community.huayno.interface import Huayno

huayno_labels = ["SHARED2", "CC2", "CC2_KEPLER", "OK2", "SHARED4"]
huayno_inttypes = [
    Huayno.inttypes.SHARED2,
    Huayno.inttypes.CC,
    Huayno.inttypes.CC_KEPLER,
    Huayno.inttypes.OK,
    Huayno.inttypes.SHARED4,
]
eta = 0.01
reference_labels = ["Hermite", "PhiGRAPE"]
all_labels = huayno_labels + reference_labels

colors = nbe._getColors()
markers = nbe._getMarkers()

# integrate the problem for a total time of T, dividing this into n_samples chunks of equal size
T = 10.0 * 2*math.pi | nbu.time
n_samples = 32.0
def defNBodyComputation():
    return nbe.NBodyComputation(
        initialConditions = mn08.planar_crisscross_initial_conditions(), analyticSolution = mn08.planar_crisscross_analytic_solution,
        dt = T / float(n_samples), tFinal = T, ndim = 2, outfName = "criss_cross"
    )

def compute():
    rw = ce.ResultsWriter() # initialize logging
    r = defNBodyComputation() # initialize the problem
    # Huayno integrators
    for (label, inttype) in zip(huayno_labels, huayno_inttypes):
        nb = Huayno(redirection='none')
        nb.initialize_code()
        nb.set_inttype_parameter(inttype)
        nb.set_timestep_parameter(eta)
        nb.set_eps2(0 | nbu.length **2)
        r.runProblemOnIntegrator(nb, label=label, printProgress = True, storeHuaynoStats = False)
    # reference integrators
    for label in reference_labels:
        nb = None
        if (label == "PhiGRAPE"):
            nb = PhiGRAPE()
            nb.set_eta(eta=eta, etas=eta)
        if (label == "Hermite"):
            nb = Hermite()
            nb.set_dt_param(eta)
        nb.set_eps2(0 | nbu.length **2)
        nb.initialize_code()
        r.runProblemOnIntegrator(nb, label = label, printProgress = True, storeHuaynoStats = False)
    r.saveResultsToDisk(rw) # save all results to disk
    del rw # stop logging

def visualize():
    # load results
    rr = ce.ResultsReader()
    r = defNBodyComputation()
    r.loadResultsFromDisk(rr)
    rr.cdto_expdir(subDir="_fig")
    # plot conserved quantities
    fig = plt.figure()
    fig.set_size_inches( (12, 8) )
    ax = fig.add_subplot(321)
    r.plotAxTotalEnergyRelativeError(ax, all_labels, logYScale=True, colors=colors, markers=markers)
    ax = fig.add_subplot(322)
    r.plotAxMomentum(ax, all_labels, logYScale=True, colors=colors, markers=markers)
    ax = fig.add_subplot(323)
    r.plotAxCOMVel(ax, all_labels, logYScale=True, colors=colors, markers=markers)
    ax = fig.add_subplot(324)
    r.plotAxAngularMomentum(ax, all_labels, logYScale=True, colors=colors, markers=markers)
    ax = fig.add_subplot(325)
    r.plotAxTrajError(ax, all_labels, logYScale=True, colors=colors, markers=markers)
    ax.legend(all_labels, bbox_to_anchor=(1.6, 1.0))
    fig.savefig("conservation-%s.pdf" % (r.outfName))
    rr.cdto_prevdir()

if __name__ == '__main__':
    # commenting compute() will run visualize() on the last results directory
    compute()
    visualize()
    print "Aloha!"
