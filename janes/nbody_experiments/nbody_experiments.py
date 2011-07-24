# -*- coding: utf-8 -*-
"""
Helper classes for studying the behaviour of N-body integrators.
Jürgen Jänes <jurgen.janes@gmail.com>

Warning: This package is mscware and needs significant refactoring and documentation
to become effortlessly usable.

A rough outline on the intended way the classes should go together follows.

    // initialize the problem
    r = NBodyComputation(
        initialConditions = <particles>, 
        analyticSolution = <function>,
        dt = 0.5 | nbody_system.time,
        tFinal = 100 | nbody_system.time,
        storeTrajectories=True
    )
    
    // run integrators on the problem
    nb = <...>
    r.runProblemOnIntegrator(nb, label="phiGrape")
    nb = <...>
    r.runProblemOnIntegrator(nb, label="Hermite")
    
    // these also return the plot objects
    f = r.plotEnergyRelativeError(labels = ["phiGrape", "Hermite"]) # returns the figure object! 
    f.savefig("...")
    r.plotImpulseRelativeError()
    r.plotTrajectories(labels = [""], plotAnalyticSolution = True)
    r.plotCoordinates(coordinates = [[0, 0], [0, 1], [0, 2]]) // plots the coordinate
    r.plotCoordinateError(coordinates = [[0, 0], [0, 1], [0, 2]]) // plots coordinate error, needs analytic solution

The 'problems' sub-package contains initial conditions as well as analytic solutions for several test problems, including:

    standard two-body problem
    Cuboctahedron orbit from [MN08]
    Planar criss-cross orbit from [MN08]

"""

import copy
import sys
import os
import math

import numpy as np
#import scipy as sp

import time

from amuse.support.io import read_set_from_file, write_set_to_file

from amuse.support.units import units
from amuse.support.units import nbody_system
import amuse.support.units.nbody_system as nbu
from amuse.support.data import core

try:
    import matplotlib
    matplotlib.use("PDF")
    import matplotlib.transforms
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d, Axes3D 
    HAS_MATPLOTLIB = True

except ImportError:
    HAS_MATPLOTLIB = False

def _getLinestyles():
    return ['-', '--', '-.', ':']

def _getColors():
    return ['k', 'y', 'b', 'm', 'c']
    #return ['0.00', '0.25', '0.50', '0.75'] #grayscale

def _getMarkers():
    from matplotlib.lines import Line2D
    """
    markers = []
    for m in Line2D.markers:
        try:
            if len(m) == 1 and m != ' ':
                markers.append(m)
        except TypeError:
            pass
    """
    # matplotlib marker symbols, re-arranged by more visible ones first
    markers = ['D', 's', 'o', 'd', '^',  'h',  '*', ',', '.', '1', 'p', '3', '2', '4', 'H', 'v', 'x', '<', '>', '|', '+', '_', ]
    return markers

def _scaleByFirst(a):
    if len(a) == 0: 
        return []
    
    b = []
    for i in range(len(a)):
        b.append((a[i] - a[0]) / a[0])
    return b

def _vprod(a,b):
    return [ 
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]
        ]

def particles_from_floats(m, x, v):
    stars = core.Stars(len(x))
    for i in range(len(stars)):
        star = stars[i]
        star.mass = m[i] | nbody_system.mass
        star.position = x[i] | nbody_system.length
        star.velocity = v[i] | ( nbody_system.length / nbody_system.time )
        star.radius = 0. | nbody_system.length # obligatory for Hermite() integrator
    return stars

def all_timestep_indices(n = 5):
    ts_i = []
    for i in range(n):
        for j in range(0, i): 
            ts_i.append( (i, j))
    return ts_i

class NBodyComputationResult:
    """
    Auxiliary class for storing the results of a computation with an N-body integrator.
    Would probably be more elegant by using the snapshot-functions implemented in the Particles-class.
    """

    def __init__(self, nb=None, storeHuaynoStats=False, storeTimesteps=[]):
        if not (nb is None):
            self.time = [ nb.get_time().value_in(nbody_system.time) ]
            self.systemStates = [ nb.particles.copy() ]
            self.kineticEnergy = [ nb.get_kinetic_energy() ]
            self.potentialEnergy = [ nb.get_potential_energy() ]
            self.totalEnergy = [ self.kineticEnergy[-1] + self.potentialEnergy[-1] ]
            self.totalEnergyRelativeError = [ 0 ]
            self.wallClockTime = [ ]

        if storeHuaynoStats:
            self.ttot = [];
            self.ktot = [];
            self.dtot = [];
            self.tstot = [];
            self.kstot = [];
            self.dstot = [];
            self.cetot = [];
            self.cetotfail= [];
        
        if len(storeTimesteps) > 0:
            self.timesteps = []
            self.timesteps_adj = []
            for i in range(len(storeTimesteps)):
                self.timesteps.append([]) 
                self.timesteps_adj.append([]) 

    def appendState(self, nb, storeHuaynoStats, wallClockTime):
        self.time.append( nb.get_time().value_in(nbody_system.time) )
        self.systemStates.append( nb.particles.copy() )
        self.kineticEnergy.append( nb.get_kinetic_energy() )
        self.potentialEnergy.append( nb.get_potential_energy() )
        #self.kineticEnergy.append( nb.particles.kinetic_energy() )
        #print "kineticEnergy: %E" % (self.kineticEnergy[-1].value_in(nbody_system.length**2 * nbody_system.time**-2 * nbody_system.mass), )
        #self.potentialEnergy.append( nb.particles.potential_energy(G=1) )
        #print "potentialEnergy: %E" % (self.potentialEnergy[-1].value_in(nbody_system.length**-1 * nbody_system.mass**2), )
        self.totalEnergy.append( self.kineticEnergy[-1] + self.potentialEnergy[-1] )
        self.totalEnergyRelativeError.append(((self.totalEnergy[-1] - self.totalEnergy[0]) / self.totalEnergy[0]).value_in(units.none))
        self.wallClockTime.append( wallClockTime )

        if storeHuaynoStats:
            stats = nb.get_evolve_statistics()
            self.ttot.append(stats[0].value_in(units.none));
            self.ktot.append(stats[1].value_in(units.none));
            self.dtot.append(stats[2].value_in(units.none));
            self.tstot.append(stats[3].value_in(units.none));
            self.kstot.append(stats[4].value_in(units.none));
            self.dstot.append(stats[5].value_in(units.none));
            self.cetot.append(stats[6].value_in(units.none));
            self.cetotfail.append(stats[7].value_in(units.none));

    def appendTimesteps(self, ts, ts_adj):
        for i in range(len(ts)):
            self.timesteps[i].append( ts[i] )
            self.timesteps_adj[i].append( ts_adj[i] )
    
    def getAllTimestepsAtTimepoint(self, timepoint):
        ts = []
        # self.timesteps is a list of timeseries with time steps
        for parti in range(len(self.timesteps)):
            ts.append(self.timesteps[parti][timepoint])
        return ts

    def getAllAdjTimestepsAtTimepoint(self, timepoint):
        ts_adj = []
        # self.timesteps is a list of timeseries with time steps
        for parti in range(len(self.timesteps_adj)):
            ts_adj.append(self.timesteps_adj[parti][timepoint])
        return ts_adj

    def getCoordinateTimeline(self, particleIndex, axisName):
        traj = []        
        for t in range(len(self.systemStates)):
            currentState = self.systemStates[t]
            particleKey = currentState.get_all_keys_in_store()[particleIndex]
            traj.append( currentState.get_value_in_store(particleKey, axisName).value_in(nbody_system.length) )
        return traj

    def _getParticlesMomentum(self, particles, axisName):
        momentum = 0 | nbody_system.mass * nbody_system.length / nbody_system.time

        for particle in particles:
            #x = particles.get_value_in_store(particle.key, axisName)
            vx = particles.get_value_in_store(particle.key, "v" + axisName)
            momentum += particle.mass * vx
            #print "%e" % (momentum.value_in(nbody_system.mass * nbody_system.length / nbody_system.time),)
            
        return momentum.value_in(nbody_system.mass * nbody_system.length / nbody_system.time)

    def getMomentumTimeline(self, axisName):
        momentumTimeline = []
        
        for t in range(len(self.systemStates)):
            momentumTimeline.append(self._getParticlesMomentum(self.systemStates[t], axisName))

        return momentumTimeline

    def _getParticlesAngularMomentum(self, particles):
        angular_momentum_unit = nbody_system.length**2 * nbody_system.time**-1 * nbody_system.mass 
        angular_momentum = [0 | angular_momentum_unit, 0 | angular_momentum_unit, 0 | angular_momentum_unit]
        for particle in particles:
            #print _vprod([particle.x, particle.y, particle.z], 
            #    [particle.mass * particle.vx, particle.mass * particle.vy, particle.mass * particle.vz])
            #print angular_momentum
            diff = _vprod([particle.x, particle.y, particle.z], 
                          [particle.mass * particle.vx, particle.mass * particle.vy, particle.mass * particle.vz])
            angular_momentum[0] += diff[0] 
            angular_momentum[1] += diff[1] 
            angular_momentum[2] += diff[2] 
        return [angular_momentum[0].value_in(angular_momentum_unit),
                angular_momentum[1].value_in(angular_momentum_unit),
                angular_momentum[2].value_in(angular_momentum_unit)]

    def getAngularMomentumTimeline(self):
        angularMomentumTimeline = []
        for t in range(len(self.systemStates)):
            angularMomentumTimeline.append(self._getParticlesAngularMomentum(self.systemStates[t]))
        return angularMomentumTimeline

    def getCOMPosTimeline(self, axisName):
        comPOSTimeline = []
        axis2index = dict(zip(['x', 'y', 'z'], range(3)))
        for t in range(len(self.systemStates)):
            com = self.systemStates[t].center_of_mass()
            comPOSTimeline.append(com[axis2index[axisName]].value_in(nbody_system.length))
        return comPOSTimeline

    def getCOMVelTimeline(self, axisName):
        comVelTimeline = []
        axis2index = dict(zip(['x', 'y', 'z'], range(3)))
        for t in range(len(self.systemStates)):
            com = self.systemStates[t].center_of_mass_velocity()
            comVelTimeline.append(com[axis2index[axisName]].value_in(nbody_system.length / nbody_system.time))
        return comVelTimeline#_scaleByFirst(comVelTimeline)

    def saveToDisk(self, label, rw):
        rw.cdto_expdir()
        for (i, systemState) in zip(range(len(self.systemStates)), self.systemStates):
            write_set_to_file(systemState, "%s-systemStates-%d.hdf5" % (label, i), "amuse")
        rw.cdto_prevdir()

        systemStatesTemp = self.systemStates
        self.systemStates = [None] * len(systemStatesTemp)
        rw.save_variable(self, label)
        self.systemStates = systemStatesTemp

    def loadFromDisk(self, label, rr):
        loaded = rr.load_variable(label)
        # http://stackoverflow.com/questions/243836/how-to-copy-all-properties-of-an-object-to-another-object-in-python
        self.__dict__ = loaded.__dict__.copy()

        rr.cdto_expdir()
        for i in range(len(self.systemStates)):
            self.systemStates[i] = read_set_from_file("%s-systemStates-%d.hdf5" % (label, i), "amuse")
            self.systemStates[i].copy_to_memory()
        rr.cdto_prevdir()

class NBodyComputation:

    def __init__(self, initialConditions, analyticSolution, dt, tFinal, ndim, storeParticles=True, outfName="nbody"):

        # carry arguments
        self.initialConditions = initialConditions.copy()
        self.analyticSolution = analyticSolution
        self.dt = dt
        self.tFinal = tFinal
        self.ndim = ndim
        self.storeParticles = storeParticles
        self.outfName = outfName
        
        # dictionary containing NBodyComputationResult's with labels as keys
        self.results = {}

    def saveResultsToDisk(self, rw):
        rw.save_variable(self.results.keys(), self.outfName)
        for label, result in self.results.items():
            result.saveToDisk("%s-%s" % (self.outfName, label), rw)

    def loadResultsFromDisk(self, rr):
        result_keys = rr.load_variable(self.outfName)
        for key in result_keys:
            self.results[key] = None
            
        for label in self.results.keys():
            self.results[ label ] = NBodyComputationResult()
            self.results[ label ].loadFromDisk("%s-%s" % (self.outfName, label), rr)

    def runProblemOnIntegrator(self, nb, label, printProgress = True, storeHuaynoStats = False, storeTimesteps = []):
        """
        Run the integrator nb on the initial conditions (defined in the constructor),
        storing the results under the key /label/.
        """
        
        # add particles via the high-level interface
        stars = self.initialConditions.copy()
        nb.particles.add_particles(stars)
        nb.commit_particles()
        from_model_to_gravity = stars.new_channel_to(nb.particles)
        from_gravity_to_model = nb.particles.new_channel_to(stars)
    
        res = NBodyComputationResult(nb, storeHuaynoStats, storeTimesteps)
        ts_ok = []
        ts_rok = []
        for (i,j) in storeTimesteps:
            ts_ok.append(nb.get_ok_timestep_ij_fw(i, j)['timestep'])
            ts_rok.append(nb.get_rok_timestep_ij_fw(i, j)['timestep'])
        res.appendTimesteps(ts_ok, ts_rok) 
        
        t = 0 | nbody_system.time
        while abs(t) < abs(self.tFinal):
            
            if printProgress:
                evolve_msg = "evolving '%s': t = %f tFinal = %f " % (label, t.value_in(nbody_system.time), self.tFinal.value_in(nbody_system.time))
                print evolve_msg
    
            ts_ok = []
            ts_rok = []
            for (i,j) in storeTimesteps:
                ts_ok.append(nb.get_ok_timestep_ij_fw(i, j)['timestep'])
                ts_rok.append(nb.get_rok_timestep_ij_fw(i, j)['timestep'])
            res.appendTimesteps(ts_ok, ts_rok) 

            t = t + self.dt
            wallClock1 = time.time()
            nb.evolve_model(t)
            wallClock2 = time.time()
            wallClockTime = wallClock2 - wallClock1
            
            #(e, px, py, px, Lx, Ly, Lz, cx, cy, cz, error) = nb.get_conserved_quantities()
            #print "%E %E %E %E %E %E %E %E %E %E %d" % (e, px, py, px, Lx, Ly, Lz, cx, cy, cz, error)

            res.appendState(nb, storeHuaynoStats, wallClockTime)

            if printProgress: print "%s time = %s sec" % (evolve_msg, wallClockTime)

        nb.stop()
        self.results[ label ] = res

    def _writeFigure(self, fig, subname, tag=None):
        fname = "%s-%s%s.pdf" % (self.outfName, subname, '' if (tag is None) else "-%s" % (tag,))
        fig.savefig(fname)
        print "NBodyComputation: wrote figure %s" % (fname,)
            
    def plotTotalEnergyRelativeError(self, labels, writeFigure=True, logYScale=False, tag=None, extAx=None, plotLegend=True):

        if extAx is None:
            # plot relative error
            fig = plt.figure()
            #fig.suptitle("XY")
            fig.set_size_inches( (12, 6) )
            ax = fig.add_subplot(111)
        else:
            ax = extAx

        ax.grid(True)

        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            if not logYScale:
                #print res.time
                #print res.totalEnergyRelativeError
                ax.plot(res.time, res.totalEnergyRelativeError, linestyle='', marker=markers[i])
            else:
                #print "rel energy log plot: (... %E, %E, %E)" % (res.totalEnergyRelativeError[-3], res.totalEnergyRelativeError[-2], res.totalEnergyRelativeError[-1])
                ax.semilogy(res.time, map(abs, res.totalEnergyRelativeError), linestyle='', marker=markers[i])

        ax.set_xlabel('time')
        ax.set_ylabel('total energy, relative error')

        if plotLegend:
            leg = ax.legend(labels, 'best')
            for t in leg.get_texts():
                t.set_fontsize('small') 

        if extAx is None:
            if (writeFigure): self._writeFigure(fig, 'energyRelError', tag)
            return fig

    def plotMomentum(self, labels, writeFigure=True, tag=None):
    
        # plot relative error
        fig = plt.figure()
        fig.suptitle("time vs total momentum on (x, y, z)-axis")
        fig.set_size_inches( (12, 7) )

        ax = fig.add_subplot(311)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            ax.plot(res.time, res.getMomentumTimeline("x"), linestyle='', marker=markers[i])
        #ax.set_xlabel('time')
        ax.set_ylabel('$p_x$')
        leg = ax.legend(labels, 'best')
        for t in leg.get_texts():
            t.set_fontsize('small') 

        ax = fig.add_subplot(312)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            ax.plot(res.time, res.getMomentumTimeline("y"), linestyle='', marker=markers[i])
        #ax.set_xlabel('time')
        ax.set_ylabel('$p_y$')
        leg = ax.legend(labels, 'best')
        for t in leg.get_texts():
            t.set_fontsize('small') 

        if (self.ndim == 3):
            ax = fig.add_subplot(313)
            ax.grid(True)
            markers = _getMarkers()
            for i in range(len(labels)):
                label = labels[i]
                res = self.results[ label ]
                ax.plot(res.time, res.getMomentumTimeline("z"), linestyle='', marker=markers[i])
            ax.set_xlabel('time /nbody time units/')
            ax.set_ylabel('$p_z$')
            leg = ax.legend(labels, 'best')
            for t in leg.get_texts():
                t.set_fontsize('small') 

        if (writeFigure): self._writeFigure(fig, 'momentum', tag)
        return fig

    def plotCOMPosition(self, labels, writeFigure=True, logYScale=False, tag=None):
    
        # plot relative error
        fig = plt.figure()
        fig.suptitle("time vs center of mass coordinate of (x, y, z)-axis")
        fig.set_size_inches( (12, 7) )

        ax = fig.add_subplot(311)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            if not logYScale:
                ax.plot(res.time, res.getCOMPosTimeline("x"), linestyle='', marker=markers[i])
            else:
                #print "rel energy log plot: (... %E, %E, %E)" % (res.totalEnergyRelativeError[-3], res.totalEnergyRelativeError[-2], res.totalEnergyRelativeError[-1])
                ax.semilogy(res.time, map(abs, res.getCOMPosTimeline("x")), linestyle='', marker=markers[i])            
        #ax.set_xlabel('time')
        ax.set_ylabel('$R_x$')
        leg = ax.legend(labels, 'best')
        for t in leg.get_texts():
            t.set_fontsize('small') 

        ax = fig.add_subplot(312)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            if not logYScale:
                ax.plot(res.time, res.getCOMPosTimeline("y"), linestyle='', marker=markers[i])
            else:
                #print "rel energy log plot: (... %E, %E, %E)" % (res.totalEnergyRelativeError[-3], res.totalEnergyRelativeError[-2], res.totalEnergyRelativeError[-1])
                ax.semilogy(res.time, map(abs, res.getCOMPosTimeline("z")), linestyle='', marker=markers[i])            
        #ax.set_xlabel('time')
        ax.set_ylabel('$R_y$')
        leg = ax.legend(labels, 'best')
        for t in leg.get_texts():
            t.set_fontsize('small') 

        if (self.ndim == 3):
            ax = fig.add_subplot(313)
            ax.grid(True)
            markers = _getMarkers()
            for i in range(len(labels)):
                label = labels[i]
                res = self.results[ label ]
                if not logYScale:
                    ax.plot(res.time, res.getCOMPosTimeline("z"), linestyle='', marker=markers[i])
                else:
                    #print "rel energy log plot: (... %E, %E, %E)" % (res.totalEnergyRelativeError[-3], res.totalEnergyRelativeError[-2], res.totalEnergyRelativeError[-1])
                    ax.semilogy(res.time, map(abs, res.getCOMPosTimeline("z")), linestyle='', marker=markers[i])            
            ax.set_xlabel('time /nbody time units/')
            ax.set_ylabel('$R_z$')
            leg = ax.legend(labels, 'best')
            for t in leg.get_texts():
                t.set_fontsize('small') 

        if (writeFigure): self._writeFigure(fig, 'com-pos', tag)
        return fig
    
    def plotCOMVelocity(self, labels, writeFigure=True, logYScale=False, tag=None):
    
        # plot relative error
        fig = plt.figure()
        fig.suptitle("time vs center of mass velocity of (x, y, z)-axis, relative error")
        fig.set_size_inches( (12, 7) )

        ax = fig.add_subplot(311)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            if not logYScale:
                ax.plot(res.time, res.getCOMVelTimeline("x"), linestyle='', marker=markers[i])
            else:
                #print "rel energy log plot: (... %E, %E, %E)" % (res.totalEnergyRelativeError[-3], res.totalEnergyRelativeError[-2], res.totalEnergyRelativeError[-1])
                ax.semilogy(res.time, map(abs, res.getCOMVelTimeline("x")), linestyle='', marker=markers[i])            
        #ax.set_xlabel('time')
        ax.set_ylabel('$R_x$')
        leg = ax.legend(labels, 'best')
        for t in leg.get_texts():
            t.set_fontsize('small') 

        ax = fig.add_subplot(312)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            if not logYScale:
                ax.plot(res.time, res.getCOMVelTimeline("y"), linestyle='', marker=markers[i])
            else:
                #print "rel energy log plot: (... %E, %E, %E)" % (res.totalEnergyRelativeError[-3], res.totalEnergyRelativeError[-2], res.totalEnergyRelativeError[-1])
                ax.semilogy(res.time, map(abs, res.getCOMVelTimeline("z")), linestyle='', marker=markers[i])            
        #ax.set_xlabel('time')
        ax.set_ylabel('$R_y$')
        leg = ax.legend(labels, 'best')
        for t in leg.get_texts():
            t.set_fontsize('small') 

        if (self.ndim == 3):
            ax = fig.add_subplot(313)
            ax.grid(True)
            markers = _getMarkers()
            for i in range(len(labels)):
                label = labels[i]
                res = self.results[ label ]
                if not logYScale:
                    ax.plot(res.time, res.getCOMVelTimeline("z"), linestyle='', marker=markers[i])
                else:
                    #print "rel energy log plot: (... %E, %E, %E)" % (res.totalEnergyRelativeError[-3], res.totalEnergyRelativeError[-2], res.totalEnergyRelativeError[-1])
                    ax.semilogy(res.time, map(abs, res.getCOMVelTimeline("z")), linestyle='', marker=markers[i])            
            ax.set_xlabel('time /nbody time units/')
            ax.set_ylabel('$R_z$')
            leg = ax.legend(labels, 'best')
            for t in leg.get_texts():
                t.set_fontsize('small') 

        if (writeFigure): self._writeFigure(fig, 'com-vel', tag)
        return fig

    def plotTrajectories(self, label, writeFigure=True):
        """
        Plot trajectories of all particles for a single integrator.
        """        
        fig = plt.figure()

        if (self.ndim == 2):
            ax = fig.add_subplot(111, aspect='equal')
    
        elif (self.ndim == 3):
            # http://stackoverflow.com/questions/3810865/need-help-with-matplotlib
            ax = Axes3D(fig)

        fig.suptitle("Trajectories for integrator: %s" % (label,))

        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        if (self.ndim >= 3): # should we ever need to plot in more than 3 dimensions...
            ax.set_zlabel('Z coordinate')

        colors = ('r', 'g', 'b', 'k')
        res = self.results[label]

        if not (self.analyticSolution is None):
            analytic_sol = self.analyticSolution( res.time )

        for i in range(len(res.systemStates[0])):
            col = colors[i]
            
            if (self.ndim == 2):
                ax.scatter(res.getCoordinateTimeline(i, 'x'), res.getCoordinateTimeline(i, 'y'), 
                           color = col, marker='x')
                if not (self.analyticSolution is None):
                    ax.plot(analytic_sol[i][0], analytic_sol[i][1], color='k')
    
            elif (self.ndim == 3):
                ax.scatter(res.getCoordinateTimeline(i, 'x'), res.getCoordinateTimeline(i, 'y'), 
                    res.getCoordinateTimeline(i, 'z'), color=col, marker='x')
                if not (self.analyticSolution is None):
                    ax.plot(analytic_sol[i][0], analytic_sol[i][1], analytic_sol[i][2], color='k')

        if (writeFigure):
            fname = "%s-traj-%s.pdf" % (self.outfName, label)
            fig.savefig(fname)
            print "plotTrajectories: wrote figure %s" % (fname,)

        return fig

    def plotHuaynoStatsCounts(self, labels, writeFigure=True, tag=None):
    
        # plot relative error
        fig = plt.figure()
        fig.suptitle("time vs formula evaluations")
        fig.set_size_inches( (12, 7) )

        ax = fig.add_subplot(311)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            ax.plot(res.time, res.ttot, linestyle='', marker=markers[i])
        #ax.set_xlabel('time')
        ax.set_ylabel('time steps')
        leg = ax.legend(labels, 'best')
        for t in leg.get_texts():
            t.set_fontsize('small') 

        ax = fig.add_subplot(312)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            ax.plot(res.time, res.ktot, linestyle='', marker=markers[i])
        #ax.set_xlabel('time')
        ax.set_ylabel('kicks')
        leg = ax.legend(labels, 'best')
        for t in leg.get_texts():
            t.set_fontsize('small') 

        ax = fig.add_subplot(313)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            ax.plot(res.time, res.dtot, linestyle='', marker=markers[i])
        ax.set_xlabel('time')
        ax.set_ylabel('drifts')
        leg = ax.legend(labels, 'best')
        for t in leg.get_texts():
            t.set_fontsize('small') 

        if (writeFigure): self._writeFigure(fig, 'huaynoCounts', tag)
        return fig
    
    def plotHuaynoStatsSteps(self, labels, writeFigure=True, tag=None):
    
        # plot relative error
        fig = plt.figure()
        fig.suptitle("time vs no of calls")
        fig.set_size_inches( (12, 7) )

        ax = fig.add_subplot(311)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            ax.plot(res.time, res.tstot, linestyle='', marker=markers[i])
        #ax.set_xlabel('time')
        ax.set_ylabel('time steps')
        leg = ax.legend(labels, 'best')
        for t in leg.get_texts():
            t.set_fontsize('small') 

        ax = fig.add_subplot(312)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            ax.plot(res.time, res.kstot, linestyle='', marker=markers[i])
        #ax.set_xlabel('time')
        ax.set_ylabel('kicks')
        leg = ax.legend(labels, 'best')
        for t in leg.get_texts():
            t.set_fontsize('small') 

        ax = fig.add_subplot(313)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            ax.plot(res.time, res.dstot, linestyle='', marker=markers[i])
        ax.set_xlabel('time')
        ax.set_ylabel('drifts')
        leg = ax.legend(labels, 'best')
        for t in leg.get_texts():
            t.set_fontsize('small') 

        if (writeFigure): self._writeFigure(fig, 'huaynoSteps', tag)
        return fig

    def plotWallClockTime(self, labels, writeFigure=True, tag=None):
    
        # plot relative error
        fig = plt.figure()
        fig.suptitle("wall clock time per time step")
        #fig.set_size_inches( (12, 7) )

        ax = fig.add_subplot(111)
        ax.grid(True)
        markers = _getMarkers()
        for i in range(len(labels)):
            label = labels[i]
            res = self.results[ label ]
            ax.plot(res.time[1:], res.wallClockTime, linestyle='', marker=markers[i])
        ax.set_xlabel('time')
        ax.set_ylabel('wall-clock time')
        leg = ax.legend(labels)
        for t in leg.get_texts():
            t.set_fontsize('small') 

        if (writeFigure): self._writeFigure(fig, 'wallClockTime', tag)
        return fig

    def plotAxQuantityTimeline(self, ax, labels, getQuantity, logYScale=False, colors=None, markers=None, refLabel=None, linestyle=None):
        if colors is None: colors = ['k'] * len(labels)
        if markers is None: markers = _getMarkers()
        if linestyle is None: linestyle = '-'
        ax.grid(True)
        for (label, marker, color) in zip(labels, markers, colors):
            res = getQuantity(label)
            x = res[0]
            y = res[1]
            if (len(x) == len(y) + 1):
                x = x[1:]
            #print len(x)
            #print len(y)
            if not (refLabel is None):
                y = np.array(y) / np.array(getQuantity(refLabel)[1])
            if not logYScale:
                ax.plot(x, y, linestyle=linestyle, color=color)
            else:
                ax.semilogy(x, y, marker=marker, linestyle=linestyle, color=color)
    
    def plotAxTotalEnergyRelativeError(self, ax, labels, **kwargs):
        """
        print "plotAxTotalEnergyRelativeError: labels=%s" % (labels,)
        for label in labels:
            print "%s %s %s" % (label, 
                                len(self.results[label].time),
                                len(self.results[label].totalEnergyRelativeError))
            print self.results[label].totalEnergyRelativeError
        """
        #print map(abs, self.results["PhiGRAPE"].totalEnergy)
        #print map(abs, self.results["Hermite"].totalEnergy)
        self.plotAxQuantityTimeline(ax, labels, 
            lambda label: (self.results[label].time, map(abs, self.results[label].totalEnergyRelativeError)),
            **kwargs)
        ax.set_ylabel('$|(E(t)-E(0))/E(0)|$')
    
    def plotAxMomentum(self, ax, labels, **kwargs):
        def _plotMomentumMagnitude(label):
            px0 = np.array(self.results[label].getMomentumTimeline('x'))[0]
            py0 = np.array(self.results[label].getMomentumTimeline('y'))[0]
            pz0 = np.array(self.results[label].getMomentumTimeline('z'))[0]
            px = np.array(self.results[label].getMomentumTimeline('x'))
            py = np.array(self.results[label].getMomentumTimeline('y'))
            pz = np.array(self.results[label].getMomentumTimeline('z'))
            p = (px-px0)**2 + (py-py0)**2 + (pz-pz0)**2
            return (self.results[label].time, map(math.sqrt, p))
        self.plotAxQuantityTimeline(ax, labels, _plotMomentumMagnitude, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('$|\mathbf{p}(t)-\mathbf{p}(0)|$')
    
    def plotAxCOMPos(self, ax, labels, **kwargs):
        def _plotComPosMagnitude(label):
            x0 = np.array(self.results[label].getCOMPosTimeline('x'))[0]
            y0 = np.array(self.results[label].getCOMPosTimeline('y'))[0]
            z0 = np.array(self.results[label].getCOMPosTimeline('z'))[0]
            x = np.array(self.results[label].getCOMPosTimeline('x'))
            y = np.array(self.results[label].getCOMPosTimeline('y'))
            z = np.array(self.results[label].getCOMPosTimeline('z'))
            r = (x-x0)**2 + (y-y0)**2 + (z-z0)**2
            return (self.results[label].time, map(math.sqrt, r))
        self.plotAxQuantityTimeline(ax, labels, _plotComPosMagnitude, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('$|\mathbf{x}_{cm}(t)-\mathbf{x}_{cm}(0)|$')

    def plotAxCOMVel(self, ax, labels, **kwargs):
        def _plotComVelMagnitude(label):
            x0 = np.array(self.results[label].getCOMVelTimeline('x'))[0]
            y0 = np.array(self.results[label].getCOMVelTimeline('y'))[0]
            z0 = np.array(self.results[label].getCOMVelTimeline('z'))[0]
            x = np.array(self.results[label].getCOMVelTimeline('x'))
            y = np.array(self.results[label].getCOMVelTimeline('y'))
            z = np.array(self.results[label].getCOMVelTimeline('z'))
            r = (x-x0)**2 + (y-y0)**2 + (z-z0)**2
            return (self.results[label].time, map(math.sqrt, r))
        self.plotAxQuantityTimeline(ax, labels, _plotComVelMagnitude, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('$|\mathbf{v}_{cm}(t)-\mathbf{v}_{cm}(0)|$')
    
    def plotAxAngularMomentum(self, ax, labels, **kwargs):
        def _plotAngularMomentum(label):
            angular_momentum = self.results[label].getAngularMomentumTimeline()
            angular_momentum_mag = []
            l0 = angular_momentum[0]
            for l in angular_momentum:
                angular_momentum_mag.append(math.sqrt((l0[0]-l[0])**2 + (l0[1]-l[1])**2 + (l0[2]-l[2])**2))
            return (self.results[label].time, angular_momentum_mag)
        self.plotAxQuantityTimeline(ax, labels, _plotAngularMomentum, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('$|\mathbf{L}(t)-\mathbf{L}(0)|$')

    def plotAxTrajError(self, ax, labels, **kwargs):
        def _plotAxTrajError(label):
            res = self.results[label]
            nparticles = len( res.systemStates[0] )
            # scratch: check analytic solution drift
            #for i in range(len(res.time)):
            #    res.time[i] -= 1.00            
            analytic_sol = self.analyticSolution( res.time )
            trajError = len(res.time) * [0]
            for i in range(nparticles):
                xt = res.getCoordinateTimeline(i, 'x')
                yt = res.getCoordinateTimeline(i, 'y')
                zt = res.getCoordinateTimeline(i, 'z')
                for j in range(len(trajError)):
                    x = xt[j]
                    y = yt[j]
                    z = zt[j]
                    xa = analytic_sol[i][0][j]
                    ya = analytic_sol[i][1][j] 
                    za = analytic_sol[i][2][j]
                    err = math.sqrt((x-xa)**2 + (y-ya)**2 + (z-za)**2)
                    trajError[j] = max(trajError[j], err)                
            return (res.time, trajError)

        self.plotAxQuantityTimeline(ax, labels, _plotAxTrajError, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('$e_{\mathbf{x}}(t)$')

    def plotAxWallClockTime(self, ax, labels, **kwargs):        
        def _plotAxWallClockTime(label):
            return (self.results[label].time, self.results[label].wallClockTime)
        self.plotAxQuantityTimeline(ax, labels, _plotAxWallClockTime, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('calculation time')

    def plotAxTotalTS(self, ax, labels, **kwargs):    
        def _plotAxTotalTS(label):
            return (self.results[label].time, self.results[label].ttot)
        self.plotAxQuantityTimeline(ax, labels, _plotAxTotalTS, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('time step evaluations')

    def plotAxTotalKicks(self, ax, labels, **kwargs):
        def _plotAxTotalKicks(label):
            return (self.results[label].time, self.results[label].ktot)
        self.plotAxQuantityTimeline(ax, labels, _plotAxTotalKicks, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('kick evaluations')
        
    def plotAxTotalDrifts(self, ax, labels, **kwargs):
        def _plotAxTotalDrifts(label):
            return (self.results[label].time, self.results[label].dtot)
        self.plotAxQuantityTimeline(ax, labels, _plotAxTotalDrifts, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('drift evaluations')

    def plotAxTotalTwoBodyCalls(self, ax, labels, **kwargs):
        def _plotAxTotalDrifts(label):
            return (self.results[label].time, self.results[label].cetot)
        self.plotAxQuantityTimeline(ax, labels, _plotAxTotalDrifts, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('two body solver calls')

    def plotAxTotalTwoBodyCallsFailed(self, ax, labels, **kwargs):
        def _plotAxTotalDrifts(label):
            return (self.results[label].time, self.results[label].cetotfail)
        self.plotAxQuantityTimeline(ax, labels, _plotAxTotalDrifts, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('two body solver calls')

    def plotAxCoreDensity(self, ax, labels, **kwargs):
        def _plotAxCoreDensity(label):
            coreDensity = []
            for i in range(len(self.results[label].systemStates)):
                cd = self.results[label].systemStates[i].densitycentre_coreradius_coredens()[2].value_in(
                    nbu.mass / (nbu.length**3))
                coreDensity.append(cd)
            return (self.results[label].time, coreDensity)
        self.plotAxQuantityTimeline(ax, labels, _plotAxCoreDensity, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('core density')

    def plotAxCoreRadius(self, ax, labels, **kwargs):
        def _plotAxCoreRadius(label):
            coreRadius = []
            for i in range(len(self.results[label].systemStates)):
                cr = self.results[label].systemStates[i].densitycentre_coreradius_coredens()[1].value_in(nbu.length)
                coreRadius.append(cr)
            return (self.results[label].time, coreRadius)
        self.plotAxQuantityTimeline(ax, labels, _plotAxCoreRadius, **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('core radius')

    def plotAxPlummerRadii(self, ax, labels, **kwargs):
        def _plotAxPlummerRadius(label, mf):
            radii = []
            for i in range(len(self.results[label].systemStates)):
                r = self.results[label].systemStates[i].LagrangianRadii(mf=[mf])[0].value_in(nbu.length)
                radii.append(r)
            return (self.results[label].time, radii)
        mass_frac = [0.9, 0.5, 0.1, 0.01]
        for mf in mass_frac:
            self.plotAxQuantityTimeline(ax, labels, lambda x: _plotAxPlummerRadius(x, mf), **kwargs)
        ax.set_xlabel('')
        ax.set_ylabel('core radius')

    def plotAxTrajectories(self, ax, label, writeFigure=True):
        """
        Plot trajectories of all particles for a single integrator.
        """
        
        ax.set_xlabel('X coordinate')
        ax.set_ylabel('Y coordinate')
        if (self.ndim >= 3): # should we ever need to plot in more than 3 dimensions...
            ax.set_zlabel('Z coordinate')

        colors = ('r', 'g', 'b', 'k')
        res = self.results[label]

        for i in range(len(res.systemStates[0])):
            ax.plot([res.getCoordinateTimeline(i, 'x')[0]], [res.getCoordinateTimeline(i, 'y')[0]], 
                [res.getCoordinateTimeline(i, 'z')[0]], color='g', marker='o', linestyle='')

        # plot random trajectories
        for i in [0]:
            ax.plot([res.getCoordinateTimeline(i, 'x')[0]], [res.getCoordinateTimeline(i, 'y')[0]], 
                [res.getCoordinateTimeline(i, 'z')[0]], color='r', marker='o', linestyle='')
            ax.plot(res.getCoordinateTimeline(i, 'x'), res.getCoordinateTimeline(i, 'y'), 
                res.getCoordinateTimeline(i, 'z'), color='r', marker='', linestyle='-')

        return ax

    def plot_distance_ij_planets(self, ax = None, label=None, ndim=3, indices=None, writeFigure=False, v=False, marker='', **kwargs):
        """
        plots distances of particles (i, j) = indices
        """
        if ax is None:
            fig = plt.figure()
            fig.set_size_inches( (12, 6) )    
            ax = fig.add_subplot(111)
        res = self.results[label]
        #indices = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
        #labels = ["[0, 1]", "[0, 2]", "[0, 3]", "[1, 2]", "[1, 3]", "[2, 3]"]
        n_systems = len(res.systemStates[0]) / 2
        star_indices = range(n_systems)
        if indices is None:
            indices = []
            for i in star_indices:
                indices.append([i, i + n_systems])
        for (i, j) in indices:
            xi = np.array(res.getCoordinateTimeline(i, 'x'))
            yi = np.array(res.getCoordinateTimeline(i, 'y'))
            zi = np.array(res.getCoordinateTimeline(i, 'z'))
            xj = np.array(res.getCoordinateTimeline(j, 'x'))
            yj = np.array(res.getCoordinateTimeline(j, 'y'))
            zj = np.array(res.getCoordinateTimeline(j, 'z'))
            distij = np.sqrt((xj-xi)**2 + (yj-yi)**2 + (zj-zi)**2)
            distij = distij / distij[0]
            if v: print res.time
            if v: print distij
            ax.plot(res.time, distij, linestyle='-', marker=marker, **kwargs)
            ax.set_title("%s" % (label,))
        if (writeFigure): self._writeFigure(fig, 'distij', tag=label)
        return ax
    