import math
import optparse
import sys
import time
import numpy as np
from amuse.support.units import units
import amuse.support.units.nbody_system as nbu
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.ph4.interface import ph4
import nbody_experiments.nbody_experiments as nbe
import nbody_experiments.problems.two_body as tbp
import nbody_experiments.problems.plummer_with_planets as plummer
import nbody_experiments.problems.mn08_problems as mn08
import nbody_experiments.problems.pythagorean_test as pt
import huayno2.interface as huayno2
import ceutils.ceutils as ce
import matplotlib.pyplot as plt

hlabels2 = ["CC2_TWOBODY", "CC2", "SF2"]
inttypes2 = [
    huayno2.inttypes.CC2_TWOBODY,
    huayno2.inttypes.HOLD_DKD_CC2,
    huayno2.inttypes.HOLD_DKD,
]
hlabels4 = ["SHARED4", "CC4", "OK4"]
inttypes4 = [
    huayno2.inttypes.SHARED4,
    huayno2.inttypes.HOLD_DKD_CC4,
    huayno2.inttypes.OK4,
]
hlabels = hlabels2 #+ hlabels4
inttypes = inttypes2 #+ inttypes4
hmarks2 = [' ', ' ', ' ', ' ']
hmarks4 = ['D', 'o', 'd']
colors2 = ['k', 'r', 'g', 'b']
colors4 = ['r', 'r', 'r']
labels_all = hlabels
eta = 0.01
#indices = range(4)
indices = range(1)

def defNBodyComputation(index = 0):
    np.random.seed(index)
    n = 100 # 1000
    T = 50.0 | nbu.time #400.0 | nbu.time
    n_samples = 64.0
    return nbe.NBodyComputation(
        initialConditions = plummer.plummer_initial_conditions(n=n), analyticSolution = None,
        dt = T / n_samples, tFinal = T, ndim = 3, outfName = "plummer_cons%02d" % (index,), storeParticles = True
    )

def compute():
    rw = ce.ResultsWriter() # initialize logging
    # loop goes here!
    for index in indices:
        r = defNBodyComputation(index) # initialize the problem
        for (label, inttype) in zip(hlabels, inttypes): # Huayno integrators
            hc = huayno2.Huayno(redirection='none')
            hc.initialize_code()
            hc.set_inttype_parameter(inttype)
            hc.set_timestep_parameter(eta)
            hc.set_eps2(0.0) # since we are using close encounter support
            r.runProblemOnIntegrator(hc, label=label, printProgress = True, storeHuaynoStats = True)
        #r.runProblemOnIntegrator(PhiGRAPE(), label="PhiGRAPE", printProgress = False)
        #r.runProblemOnIntegrator(Hermite(), label="Hermite", printProgress = False)
        r.saveResultsToDisk(rw) # save all results to disk
    del rw # stop logging

def visualize():
    # load results
    rr = ce.ResultsReader()
    """
    fig_all = plt.figure()
    fig_all.set_size_inches( (12, 8) )
    ax1 = fig_all.add_subplot(221)
    ax2 = fig_all.add_subplot(222)
    ax3 = fig_all.add_subplot(223)
    ax4 = fig_all.add_subplot(224)
    """
    for index in indices:
        r = defNBodyComputation(index)
        r.loadResultsFromDisk(rr)
        rr.cdto_expdir(subDir="_fig")

        """
        # plot conservation (current index)
        fig = plt.figure()
        fig.set_size_inches( (12, 8) )
        ax = fig.add_subplot(221)
        r.plotAxTotalEnergyRelativeError(ax, hlabels2, logYScale=True, colors=colors2, markers=hmarks2)
        ax.set_xlabel('time /$N$-body units/')    
        ax = fig.add_subplot(222)
        r.plotAxMomentum(ax, hlabels2, logYScale=False, colors=colors2, markers=hmarks2)
        ax.set_xlabel('time /$N$-body units/')    
        ax = fig.add_subplot(223)
        r.plotAxCOMVel(ax, hlabels2, logYScale=False, colors=colors2, markers=hmarks2)
        ax.set_xlabel('time /$N$-body units/')    
        ax = fig.add_subplot(224)
        r.plotAxAngularMomentum(ax, hlabels2, logYScale=True, colors=colors2, markers=hmarks2)
        ax.set_xlabel('time /$N$-body units/')    
        ax.legend(hlabels2, "upper left")
        fig.savefig("conservation-%s.pdf" % (r.outfName))

        # plot dynamics
        fig = plt.figure()
        fig.set_size_inches( (12, 4) )
        ax = fig.add_subplot(121)
        r.plotAxPlummerRadii(ax, hlabels2, logYScale=True, colors=colors2, markers=hmarks2, linestyle = ':')
        r.plotAxCoreRadius(ax, hlabels2, logYScale=True, colors=colors2, markers=hmarks2)
        ax.set_xlabel('time /$N$-body units/')    
        ax.set_ylabel('{Lagrangian, core} radius')
        ax = fig.add_subplot(122)
        r.plotAxCoreDensity(ax, hlabels2, logYScale=True, colors=colors2, markers=hmarks2)
        ax.set_xlabel('time /$N$-body units/')    
        ax.set_ylabel('core density')
        ax.legend(hlabels2, 'upper left')
        fig.savefig("dynamics-%s.pdf" % (r.outfName))
        """

        # plot performance indicators (current index)
        fig = plt.figure()
        fig.set_size_inches( (12, 12) )
        ax = fig.add_subplot(321)
        #ax.set_ylim([0, 1.1])
        r.plotAxWallClockTime(ax, hlabels2, logYScale=False, colors=colors2, markers=hmarks2)#, refLabel="SF2")
        ax.set_xlabel('time /$N$-body units/')    
        ax.set_ylabel('computation time (normalized)')
        ax = fig.add_subplot(322)
        ax.set_ylim([0, 1.1])
        r.plotAxTotalTS(ax, hlabels2, logYScale=False, colors=colors2, markers=hmarks2, refLabel="SF2")
        ax.set_xlabel('time /$N$-body units/')    
        ax.set_ylabel('time step evaluations (normalized)')
        ax = fig.add_subplot(323)
        #ax.set_ylim([0, 1.1])
        r.plotAxTotalKicks(ax, hlabels2, logYScale=False, colors=colors2, markers=hmarks2, refLabel="SF2")
        ax.set_xlabel('time /$N$-body units/')    
        ax.set_ylabel('kick evaluations (normalized)')
        ax = fig.add_subplot(324)
        r.plotAxTotalDrifts(ax, hlabels2, logYScale=True, colors=colors2, markers=hmarks2, refLabel="SF2")
        ax.set_xlabel('time /$N$-body units/')    
        ax.set_ylabel('drift evaluations (normalized)')
        ax.legend(hlabels2, "lower left")
        ax = fig.add_subplot(325)
        r.plotAxTotalTwoBodyCalls(ax, hlabels2, logYScale=False, colors=colors2, markers=hmarks2)
        ax.set_xlabel('time /$N$-body units/')    
        ax.set_ylabel('Kepler solver calls')
        ax = fig.add_subplot(326)
        r.plotAxTotalTwoBodyCallsFailed(ax, hlabels2, logYScale=False, colors=colors2, markers=hmarks2)
        fig.savefig("performance-%s.pdf" % (r.outfName))
        ax.set_xlabel('time /$N$-body units/')    
        ax.set_ylabel('Kepler solver failures')
        rr.cdto_prevdir()

        """
        # plot conservation (all)
        r.plotAxTotalEnergyRelativeError(ax1, hlabels2, logYScale=True, colors=colors2, markers=hmarks2)
        ax1.set_xlabel('time /$N$-body units/')    
        r.plotAxMomentum(ax2, hlabels2, logYScale=True, colors=colors2, markers=hmarks2)
        ax2.set_xlabel('time /$N$-body units/')    
        r.plotAxCOMPos(ax3, hlabels2, logYScale=True, colors=colors2, markers=hmarks2)
        #r.plotAxVelPos(ax3, hlabels2, logYScale=True, colors=colors2, markers=hmarks2)
        ax3.set_xlabel('time /$N$-body units/')    
        r.plotAxAngularMomentum(ax4, hlabels2, logYScale=True, colors=colors2, markers=hmarks2)
        ax4.set_xlabel('time /$N$-body units/')  
        """  

    #ax.legend(labels_all, bbox_to_anchor=(2.1, 1.05))
    """
    rr.cdto_expdir(subDir="_fig")
    fig_all.savefig("conservation-plummer-all.pdf");
    rr.cdto_prevdir()
    """
if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-c", action="store_true", dest="compute_flag", default=False)
    parser.add_option("-v", action="store_true", dest="visualize_flag", default=False)
    (options, args) = parser.parse_args()
    # force default behavior to compute and visualize
    if not (options.compute_flag) and not (options.visualize_flag):
        options.compute_flag = True
        options.visualize_flag = True
    #if options.compute_flag: compute()
    if options.visualize_flag: visualize()
    print "Fin."
