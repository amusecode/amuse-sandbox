# NOTE
# You need a snapshot and plot directory:
# mkdir snapshot
# mkdir plot
# Author: Guillermo Kardolus
# Email:  kardolus@strw.leidenuniv.nl

import sys
import unittest
import numpy 
import random
import collections
import os
import re

try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.units import nbody_system
from amuse.support.units import units

from amuse.legacy.hermite0.interface import Hermite
from amuse.legacy.bhtree.interface import BHTree
from amuse.legacy.phiGRAPE.interface import PhiGRAPEInterface, PhiGRAPE
from amuse.legacy.sse.interface import SSE
from amuse.legacy.support.core import is_mpd_running

from amuse.support.io import store
from os import system

from amuse.ext.plummer import MakePlummerModel
from amuse.ext.salpeter import SalpeterIMF

bold = "\033[1m"
reset = "\033[0;0m"

system("echo \" \" > diagnostics.dat")
system("echo \" \" > mass.dat")
system("echo \" \" > lowMass.dat")
system("echo \" \" > mediumMass.dat")
system("echo \" \" > highMass.dat")
system("echo \" \" > lagrange.dat")
system("echo \" \" > simulation.log")
system("echo \" \" > binaries.dat")
system("echo \" \" > unbound.dat")
system("rm -f snapshot/* 2&>1")
system("rm -f plot/* 2&>1")
system("echo \" \" > snap.gpl")
file = open('diagnostics.dat', 'a')
massFile = open('mass.dat', 'a')
lowMassFile = open('lowMass.dat', 'a')
mediumMassFile = open('mediumMass.dat', 'a')
highMassFile = open('highMass.dat', 'a')
lagrangeFile = open('lagrange.dat', 'a')
logFile = open('simulation.log', 'a')
binaryFile = open('binaries.dat', 'a')
unboundFile = open('unbound.dat', 'a')

massFile.write('lowest_mass(MSun) highest_mass(Msun) middle_mass(MSun) NStars\n')
lowMassFile.write('lowest_mass(MSun) highest_mass(Msun) middle_mass(MSun) NStars\n')
mediumMassFile.write('lowest_mass(MSun) highest_mass(Msun) middle_mass(MSun) NStars\n')
highMassFile.write('lowest_mass(MSun) highest_mass(Msun) middle_mass(MSun) NStars\n')
file.write('time(Myr) total_mass(MSun) mass_loss(MSun) mass_loss(%) total_mass_loss(MSun) total_mass_loss(%) total_energy(J) total_energy_error(%) gravity_error(%) evolution_error(%) total_gravity_error(%) total_evolution_error(%) evo_error/gra_error kinetic/potential')
file.write('\n')
lagrangeFile.write('time(Myr) 0.5(%) 1(%) 5(%) 10(%) 25(%) 50(%) 75(%) 90(%)\n')

class SalpeterIMF(object):
    def __init__(self, mass_min = 0.1 | units.MSun, mass_max = 125 | units.MSun, alpha = -2.35):
        self.mass_min = mass_min.as_quantity_in(units.MSun)
        self.mass_max = mass_max.as_quantity_in(units.MSun)
        self.alpha = alpha
        self.random = random.Random()
        self.random.seed()
    
    def mass_mean(self):
        alpha1 = self.alpha + 1
        alpha2 = self.alpha + 2
        l1 = pow(self.mass_min.value_in(units.MSun), alpha1)
        l2 = pow(self.mass_min.value_in(units.MSun), alpha2)
        u1 = pow(self.mass_max.value_in(units.MSun), alpha1)
        u2 = pow(self.mass_max.value_in(units.MSun), alpha2)
        return ((u2 - l2) * alpha1) / ((u1 - l1) * alpha2) | units.MSun
        
    def mass(self, random_number):
        alpha1 = self.alpha + 1
        factor = (pow(self.mass_max.value_in(units.MSun) / self.mass_min.value_in(units.MSun) , alpha1) - 1.0)
        return self.mass_min.value_in(units.MSun) * (pow(1 + (factor * random_number), 1.0 / alpha1)) | units.MSun
        
    def next_mass(self):
        return self.mass(self.random.random())
        
    def next_set(self, number_of_stars):
        set_of_masses = numpy.zeros(number_of_stars)
        total_mass = 0.0 | units.MSun
        for i in range(number_of_stars):
           mass = self.next_mass()
           set_of_masses[i] = mass.value_in(units.MSun)
           total_mass += mass
        return (total_mass, units.MSun.new_quantity(set_of_masses))
        
class SalpeterIMFTests(unittest.TestCase):
    def test1(self):
        instance = SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual(instance.mass_mean().value_in(units.MSun), 0.351, 3)

    def test2(self):
        instance = SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual(instance.mass(1.0).value_in(units.MSun), 100, 3)
        self.assertAlmostEqual(instance.mass(0).value_in(units.MSun), 0.1, 3)
       
    def test3(self):
        instance = SalpeterIMF(0.1 | units.MSun, 100 | units.MSun, alpha = -2.35)
        n = 10000
        total_mass, set_of_masses = instance.next_set(10000)
        mean = total_mass.value_in(units.MSun) / float(n)
        exact_mean = instance.mass_mean().value_in(units.MSun)
        self.assertTrue(abs(mean - exact_mean) < 0.1)
        
    def test4(self):
        instance = SalpeterIMF(0.1 | units.MSun, 125 | units.MSun, alpha = -2.35)
        self.assertAlmostEqual( 1.0 / instance.mass_mean().value_in(units.MSun), 2.8253, 4)
       

def move_particles_to_center_of_mass(particles): 
    center_of_mass = particles.center_of_mass()
    center_of_mass_velocity = particles.center_of_mass_velocity()
    
    particles.position = particles.position - center_of_mass
    particles.velocity = particles.velocity - center_of_mass_velocity     
   
     
def plot_particles(particles, name_of_the_figure):
    
    if HAS_MATPLOTLIB:
        print "plotting the data"
        
        
        figure = pyplot.figure(figsize = (40, 40))
        plots = map(lambda x : figure.add_subplot(4,4, x+1), range(4*4))
        
        index = 0
        for data in particles.history:
            if index > 15:
                break
            
            x_values = data.x.value_in(units.parsec)
            y_values = data.y.value_in(units.parsec)
            mass_values = data.mass.value_in(units.MSun)
            
            sizes = mass_values * 10.0
            plots[index].scatter(x_values,y_values, marker='o', s=sizes)
            index += 1
            
        for plot in plots:
            plot.set_xlim(-2.0, 2.0)
            plot.set_ylim(-2.0, 2.0)

        figure.savefig(name_of_the_figure)
        if False:
            from matplotlib import axes3d
            figure = pyplot.figure()
            axes_3d = axes3d.Axes3D(figure)
            positions = particles.get_values_of_attribute('position')
            xs = numpy.array([position.x.value_in(units.lightyear) for position in positions])
            ys = numpy.array([position.y.value_in(units.lightyear) for position in positions])
            zs = numpy.array([position.z.value_in(units.lightyear) for position in positions])
            #print xs, yz, zs
            
            plot = axes_3d.scatter(xs, ys, zs)
            axes_3d.set_xlim(-10.0,10.0)
            axes_3d.set_ylim(-10.0,10.0)  
            axes_3d.set_zlim(-10.0,10.0)        
            
            figure.savefig("3d_"+name_of_the_figure)

def print_log(time, gravity, particles, total_energy_at_t0, total_energy_at_this_time):
    print "Evolved model to t    = " + str(time)
    print total_energy_at_t0, total_energy_at_this_time, (total_energy_at_this_time - total_energy_at_t0) / total_energy_at_t0
    print "KE:" , particles.kinetic_energy().as_quantity_in(units.J)
    print "PE:" , particles.potential_energy(gravity.parameters.epsilon_squared)
    print  "center of mass:" , particles.center_of_mass()
    print  "center of mass velocity:" , particles.center_of_mass_velocity()

    
def simulate_small_cluster(number_of_stars, end_time = 40 | units.Myr, name_of_the_figure = "test-2.svg"):
    random.seed()
    
    initial_mass_function = SalpeterIMF()
    total_mass, salpeter_masses = initial_mass_function.next_set(number_of_stars)
    
    convert_nbody = nbody_system.nbody_to_si(total_mass, 1.0 | units.parsec)
    convert_nbody.set_as_default()
    
    particles = MakePlummerModel(number_of_stars, convert_nbody).result;

    gravity = PhiGRAPE(convert_nbody,mode="gpu")
    gravity.initialize_code()
    gravity.parameters.timestep_parameter=0.01
    gravity.parameters.initial_timestep_parameter=0.01
    gravity.parameters.epsilon_squared = 0.000001 | units.parsec ** 2
     
    stellar_evolution = SSE()
    stellar_evolution.initialize_module_with_default_parameters() 
    
    particles.radius = 0.0 | units.RSun
    
    gravity.particles.add_particles(particles)
    gravity.initialize_particles(0.0)
    
    from_model_to_gravity = particles.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(particles)
   
    time = 0.0 | units.Myr    
    particles.savepoint(time)

    move_particles_to_center_of_mass(particles)

    file.write(str(time.value_in(units.Myr)))
    file.write(' ')

    #Count number of binaries
    stack = []
    nBinaries = 0
    nPotentialBinaries = 0
    nHardBinaries = 0
    averageKineticEnergy = 0.0
    for stars in particles:
      stackSize = len(stack)
      distanceSquared = pow(stars.x.value_in(units.m), 2.0) + pow(stars.y.value_in(units.m), 2.0) + pow(stars.z.value_in(units.m), 2.0)
      xPosition = stars.x.value_in(units.m)
      yPosition = stars.y.value_in(units.m)
      zPosition = stars.z.value_in(units.m)
      xVelocity = stars.vx.value_in(units.ms)
      yVelocity = stars.vy.value_in(units.ms)
      zVelocity = stars.vz.value_in(units.ms)
      lastMass = stars.mass.value_in(units.kg)      

      stack.append([stackSize])
      stack[stackSize].append(distanceSquared)  # x1 
      stack[stackSize].append(xPosition)        # x2
      stack[stackSize].append(yPosition)        # x3
      stack[stackSize].append(zPosition)        # x4
      stack[stackSize].append(lastMass)         # x5
      stack[stackSize].append(xVelocity)        # x6
      stack[stackSize].append(yVelocity)        # x7
      stack[stackSize].append(zVelocity)        # x8

      velocitySquared = pow(stars.vx.value_in(units.ms), 2.0) + pow(stars.vy.value_in(units.ms), 2.0) + pow(stars.vz.value_in(units.ms), 2.0)
      averageKineticEnergy += 0.5 * velocitySquared * stars.mass.value_in(units.kg)

    averageKineticEnergy = averageKineticEnergy / (stackSize + 1)

    #Sort stack by distance from center
    stack.sort(key=lambda x: x[1])
   
    #Look at the $nCompare closest "neighbors" (because not sure if they are actually neighbors, but if they are a binary, the odds are in favor of them being neighbors in this definition)
    nCompare = 2

    for x in range (0, stackSize):
      if x + nCompare <= stackSize:
        # While there are stars to compare with, do the computations
        for nNeighbors in range (1, nCompare + 1):
          # 1. Look at the Energy of the pairs
          y = x + nNeighbors
          xSeparationSquared = pow(stack[x][2] - stack[y][2], 2.0)
          ySeparationSquared = pow(stack[x][3] - stack[y][3], 2.0)
          zSeparationSquared = pow(stack[x][4] - stack[y][4], 2.0)
          xDeltaVelocitySquared = pow(stack[x][6] - stack[y][6], 2.0)
          yDeltaVelocitySquared = pow(stack[x][7] - stack[y][7], 2.0)
          zDeltaVelocitySquared = pow(stack[x][8] - stack[y][8], 2.0)
          separation = pow(xSeparationSquared + ySeparationSquared + zSeparationSquared, 0.5)
          deltaVelocitySquared = xDeltaVelocitySquared + yDeltaVelocitySquared + zDeltaVelocitySquared
          reducedMass = stack[x][5] * stack[y][5] / (stack[x][5] + stack[y][5]) 
          kineticEnergy = 0.5 * reducedMass * deltaVelocitySquared
          potentialEnergy = - 6.673e-11 * stack[x][5] * stack[y][5] / separation
          totalEnergy = kineticEnergy + potentialEnergy
         
          if totalEnergy < 0.0:
            # Compute Semi Major Axis (a)
            # Specific orbital energy = E_tot / mu, a = - G (m1 + m2) / 2 epsilon 
            orbitalEnergy = totalEnergy / reducedMass
            semiMajorAxis = - 6.673e-11 * (stack[x][5] + stack[y][5]) / (2 * orbitalEnergy)
            nPotentialBinaries += 1

            # 2. Look only at neighbors with a < 0.1 pc (Kroupa et al - "The formation of very wide binaries")
            if semiMajorAxis < 0.1 * 3.08568025e16:
              # 3. Look at the eccentricity
              # Calculate specific angular momentum squared (L / mu)^2 = (r x v)^2
              angularMomentum = pow((stack[x][2] - stack[y][2]) * (stack[x][7] - stack[y][7]) - (stack[x][3] - stack[y][3]) * (stack[x][6] - stack[y][6]), 2.0)
              angularMomentum += pow((stack[x][3] - stack[y][3]) * (stack[x][8] - stack[y][8]) - (stack[x][4] - stack[y][4]) * (stack[x][7] - stack[y][7]), 2.0)
              angularMomentum += pow((stack[x][4] - stack[y][4]) * (stack[x][6] - stack[y][6]) - (stack[x][2] - stack[y][2]) * (stack[x][8] - stack[y][8]), 2.0)
              eccentricity = pow(1.0 - angularMomentum / (6.673e-11 * (stack[x][5] + stack[y][5]) * semiMajorAxis), 0.5)
              if eccentricity <= 1.0 and eccentricity >= 0.0:
                nBinaries += 1
                if - potentialEnergy > averageKineticEnergy:
                  nHardBinaries +=1
    
    nSoftBinaries = nBinaries - nHardBinaries
    print "N Possible Binaries   =", nPotentialBinaries
    print "Number of Binaries    =", nBinaries
    print "N Hard Binaries       =", nHardBinaries
    print "N Soft Binaries       =", nSoftBinaries
    binaryFile.write('0.0')
    binaryFile.write(' ')
    binaryFile.write(str(nPotentialBinaries))
    binaryFile.write(' ')
    binaryFile.write(str(nBinaries))
    binaryFile.write(' ')
    binaryFile.write(str(nHardBinaries))
    binaryFile.write(' ')
    binaryFile.write(str(nSoftBinaries))
    binaryFile.write('\n')

    #Compute half-mass radius
    lagrangeMass = 0.0
    stack = []

    for stars in particles:
      stackSize = len(stack)

      distanceSquared = pow(stars.x.value_in(units.parsec), 2.0) + pow(stars.y.value_in(units.parsec), 2.0) + pow(stars.z.value_in(units.parsec), 2.0)
      lagrangeMass += stars.mass.value_in(units.MSun)
      lastMass = stars.mass.value_in(units.MSun)
      lastMassSI = stars.mass.value_in(units.kg)

      # Legenda:
      # stack[x][0] = ID | stack[x][1] = distance^2 | stack[x][2] = mass | stack[x][3] = Etot | stack[x][4] = x-position | stack[x][5] = y-position | stack[x][6] = z-position
      stack.append([stackSize])
      stack[stackSize].append(distanceSquared)
      stack[stackSize].append(lastMass)

      velocitySquared = pow(stars.vx.value_in(units.ms), 2.0) + pow(stars.vy.value_in(units.ms), 2.0) + pow(stars.vz.value_in(units.ms), 2.0)
      Etot = particles[stackSize].potential().value_in(units.ms**2) * stars.mass.value_in(units.kg) + 0.5 * velocitySquared * stars.mass.value_in(units.kg)
      stack[stackSize].append(Etot)
      stack[stackSize].append(stars.x.value_in(units.m))
      stack[stackSize].append(stars.y.value_in(units.m))
      stack[stackSize].append(stars.z.value_in(units.m))

    stackSize = len(stack)
    stack.sort(key=lambda x: x[1])
    
    # Stars sorted by distance
    # Find half mass radius
    tmpMass = 0.0
    for x in range (0, stackSize):
      tmpMass += stack[x][2]
      if tmpMass/lagrangeMass >= 0.5:
        RHalf = stack[x][1]
        print "Half-mass radius      =", RHalf, "pc"
        print "Total mass            =", lagrangeMass, "MSun"
        print "M(1/2)                =", tmpMass, "MSun"
        print "M(Tot)/M(1/2)         =", tmpMass/lagrangeMass
        break
 
    RHalf = pow(5.0, 2.0) * RHalf
    
    # Found half mass radius
    # Find stars on bigger distance than 20 * half-mass radius
    nStarsUnbound = 0.0
    nStarsOutside = 0.0
    nBound = 0.0
    realCOMx = 0.0
    realCOMy = 0.0
    realCOMz = 0.0
    realMass = 0.0
  
    for x in range (0, stackSize):
      if stack[x][1] > RHalf:
        nStarsOutside += 1.0
        if stack[x][3] > 0.0:
          nStarsUnbound += 1.0

      # Recalculate Mass
      if stack[x][1] < RHalf or stack[x][3] < 0.0:
        realMass += stars.mass.value_in(units.MSun)

    for x in range (0, stackSize):
      if stack[x][1] < RHalf or stack[x][3] < 0.0:
        #Center of Mass: R = SUM r(i) m(i) / Mtot
        realCOMx += stack[x][4] * stack[x][2] / realMass
        realCOMy += stack[x][5] * stack[x][2] / realMass
        realCOMz += stack[x][6] * stack[x][2] / realMass
        nBound += 1.0
    
    fractionUnbound = nStarsUnbound / (stackSize + 1)
    realCOM = [realCOMx, realCOMy, realCOMz] | units.m

    print "NStars R > 5 R_0.5    =", nStarsOutside
    print "Nr of unbound stars   =", nStarsUnbound
    print "Fraction unbound      =", fractionUnbound
    print "Nr of bound stars     =", nBound
    print "Fraction bound        =", nBound / (stackSize + 1)
    print "Real Center Of Mass   =", realCOM.value_in(units.parsec)

    unboundFile.write('0.0')
    unboundFile.write(' ')
    unboundFile.write(str(nStarsUnbound))
    unboundFile.write(' ')
    unboundFile.write(str(fractionUnbound))
    unboundFile.write('\n')

    #Compute Lagrangian Radii
    stack = []
    for stars in particles:
      #shift particles to center of mass
      stars.position = stars.position - realCOM
      #stars.velocity = stars.velocity - particles.center_of_mass_velocity()
      distanceSquared = pow(stars.x.value_in(units.parsec), 2.0) + pow(stars.y.value_in(units.parsec), 2.0) + pow(stars.z.value_in(units.parsec), 2.0)
      stack.append(distanceSquared)
    stackSize = len(stack)
    stack.sort()
    lagrangeFile.write(str(time.value_in(units.Myr)))
    lagrangeFile.write(' ')
    # 0.5 %
    lagrangeFile.write(str(pow(stack[int(round(len(stack)/200.0, 0))], 0.5)))
    lagrangeFile.write(' ')
    # 1 %
    lagrangeFile.write(str(pow(stack[int(round(len(stack)/100.0, 0))], 0.5)))
    lagrangeFile.write(' ')
    # 5 %
    lagrangeFile.write(str(pow(stack[int(round(len(stack)/20.0, 0))], 0.5)))
    lagrangeFile.write(' ')
    # 10 %
    lagrangeFile.write(str(pow(stack[int(round(len(stack)/10.0, 0))], 0.5)))
    lagrangeFile.write(' ')
    # 25 %
    lagrangeFile.write(str(pow(stack[int(round(len(stack)/4.0, 0))], 0.5)))
    lagrangeFile.write(' ')
    # 50 %
    lagrangeFile.write(str(pow(stack[int(round(len(stack)/2.0, 0))], 0.5)))
    lagrangeFile.write(' ')
    # 75 %
    lagrangeFile.write(str(pow(stack[int(round(len(stack)/(100.0/75.0), 0))], 0.5)))
    lagrangeFile.write(' ')
    # 90 %
    lagrangeFile.write(str(pow(stack[int(round(len(stack)/(100.0/90.0), 0))], 0.5)))
    lagrangeFile.write(' ')

    lagrangeFile.write('\n')

    #Find mass distribution
    #First find stars with highest and lowest mass

    lowestMass = stars.mass.value_in(units.MSun)
    highestMass = lowestMass

    for x in particles:

      if x.mass.value_in(units.MSun) < lowestMass:
        lowestMass = x.mass.value_in(units.MSun)

      if x.mass.value_in(units.MSun) > highestMass:
        highestMass = x.mass.value_in(units.MSun)

    print "Highest Mass          =", highestMass, "MSun"
    print "Lowest  Mass          =", lowestMass, "MSun"

    #Compute an interval
    Ninterval = 20
    interval = highestMass/Ninterval
    intervalCount = 1
    starCount = 0

    while intervalCount <= Ninterval:
      for x in particles:
        if ((intervalCount - 1) * interval) < x.mass.value_in(units.MSun) <= (intervalCount * interval):
          starCount += 1
      massFile.write(str((intervalCount - 1) * interval))
      massFile.write(' ')
      massFile.write(str(intervalCount * interval))
      massFile.write(' ')
      massFile.write(str(((intervalCount - 1) * interval + intervalCount * interval)/2))
      massFile.write(' ')
      massFile.write(str(starCount))
      massFile.write('\n')
      starCount = 0
      intervalCount += 1

    massFile.close()

    # Low Mass stars
    Ninterval = 20
    interval = 0.05
    intervalCount = 1
    starCount = 0

    while intervalCount <= Ninterval:
      for x in particles:
        if ((intervalCount - 1) * interval) < x.mass.value_in(units.MSun) <= (intervalCount * interval):
          starCount += 1
      lowMassFile.write(str((intervalCount - 1) * interval))
      lowMassFile.write(' ')
      lowMassFile.write(str(intervalCount * interval))
      lowMassFile.write(' ')
      lowMassFile.write(str(((intervalCount - 1) * interval + intervalCount * interval)/2))
      lowMassFile.write(' ')
      lowMassFile.write(str(starCount))
      lowMassFile.write('\n')
      starCount = 0
      intervalCount += 1

    lowMassFile.close()

    # Medium Mass stars
    Ninterval = 20
    interval = 0.25 
    intervalCount = 1
    starCount = 0

    while intervalCount <= Ninterval:
      for x in particles:
        if ((intervalCount - 1) * interval) < x.mass.value_in(units.MSun) <= (intervalCount * interval):
          starCount += 1
      mediumMassFile.write(str((intervalCount - 1) * interval))
      mediumMassFile.write(' ')
      mediumMassFile.write(str(intervalCount * interval))
      mediumMassFile.write(' ')
      mediumMassFile.write(str(((intervalCount - 1) * interval + intervalCount * interval)/2))
      mediumMassFile.write(' ')
      mediumMassFile.write(str(starCount))
      mediumMassFile.write('\n')
      starCount = 0
      intervalCount += 1

    mediumMassFile.close()

    # High Mass stars
    Ninterval = 20
    interval = 0.5
    intervalCount = 1
    starCount = 0

    while intervalCount <= Ninterval:
      for x in particles:
        if ((intervalCount - 1) * interval) < x.mass.value_in(units.MSun) <= (intervalCount * interval):
          starCount += 1
      highMassFile.write(str((intervalCount - 1) * interval))
      highMassFile.write(' ')
      highMassFile.write(str(intervalCount * interval))
      highMassFile.write(' ')
      highMassFile.write(str(((intervalCount - 1) * interval + intervalCount * interval)/2))
      highMassFile.write(' ')
      highMassFile.write(str(starCount))
      highMassFile.write('\n')
      starCount = 0
      intervalCount += 1

    highMassFile.close()

    total_energy_at_t0 = gravity.kinetic_energy + gravity.potential_energy
    postEvolutionEnergy = gravity.kinetic_energy.value_in(units.J) + gravity.potential_energy.value_in(units.J)
    totalGravityError = 1
    totalEvolutionError = 1
    initialMass = 0.0
    nStars = 0.0
    for x in particles:
      initialMass += x.mass.value_in(units.MSun)
      nStars += 1.0

    file.write(str(initialMass))
    file.write(' ')
    # Mass loss
    file.write('0.0 0.0 0.0 0.0 ')
    file.write(str(postEvolutionEnergy))
    file.write(' ')
    # Energy errors
    file.write('0.0 0.0 0.0 0.0 0.0 0.0')

    isBound = gravity.kinetic_energy.value_in(units.J)/gravity.potential_energy.value_in(units.J)
    file.write(' ')
    file.write(str(isBound))
    file.write('\n')

    # Set all masses to the average mass
    averageMass = initialMass / nStars
    particles.mass = averageMass | units.MSun

    lowestMass = stars.mass.value_in(units.MSun)
    highestMass = lowestMass

    for x in particles:

      if x.mass.value_in(units.MSun) < lowestMass:
        lowestMass = x.mass.value_in(units.MSun)

      if x.mass.value_in(units.MSun) > highestMass:
        highestMass = x.mass.value_in(units.MSun)

    print "------- MASS RESET --------"
    print "Highest Mass          =", highestMass, "MSun"
    print "Lowest  Mass          =", lowestMass, "MSun"

    # Write to log
    logFile.write('Average Mass = ')
    logFile.write(str(averageMass))
    logFile.write(' MSun\n')
    logFile.write('Total Initial Mass = ')
    logFile.write(str(initialMass))
    logFile.write(' MSun\n')
    #logFile.write(str(gravity.parameters.epsilon_squared | units.parsec ** 2))
    logFile.close

    Nloop = 0.0

    while time < end_time:
        time += 0.25 | units.Myr
        gravity.evolve_model(time)
        from_gravity_to_model.copy()

        preEvolutionEnergy = gravity.kinetic_energy.value_in(units.J) + gravity.potential_energy.value_in(units.J)

        totalMass = 0.0
        for x in particles:
          totalMass += x.mass.value_in(units.MSun)

        preEvolutionMass = totalMass
        gravityError = 100 * (postEvolutionEnergy - preEvolutionEnergy)/preEvolutionEnergy
        
        print "--------------------------------------------"
        print bold + "Pre-Evolution" + reset

        print "Total Energy          =", preEvolutionEnergy, "J"
        print "Total Mass            =", totalMass, "MSun"
        print "Gravity Error         =", gravityError, "%"
        print "Total Energy Error    =", 100 * (preEvolutionEnergy - total_energy_at_t0.value_in(units.J)) / total_energy_at_t0.value_in(units.J), "%"

        #stellar_evolution.evolve_model(time)

        #from_stellar_evolution_to_model.copy_attributes(["mass", "radius"])

        particles.savepoint(time)

        #from_model_to_gravity.copy_attributes(["mass"])
        postEvolutionEnergy = gravity.kinetic_energy.value_in(units.J) + gravity.potential_energy.value_in(units.J)
        totalEnergyError = 100 * (postEvolutionEnergy - total_energy_at_t0.value_in(units.J)) / total_energy_at_t0.value_in(units.J)
        evolutionError = 100 * (postEvolutionEnergy - preEvolutionEnergy) / preEvolutionEnergy

        totalMass = 0.0
        for x in particles:
          totalMass += x.mass.value_in(units.MSun)

        print bold + "Post-Evolution" + reset
        print "Total Energy          =", postEvolutionEnergy, "J"
        print "Total Mass            =", totalMass, "MSun"
        print "Total Energy Error    =", totalEnergyError, "%"
        print bold + "Post vs Pre" + reset
        print "Evolution Error       =", evolutionError, "%"
        print "Mass Loss             =", preEvolutionMass - totalMass, "MSun"
        print "Relative Mass Loss    =", (preEvolutionMass - totalMass)/preEvolutionMass * 100, "%"
        print "Total Mass Loss       =", initialMass - totalMass, "MSun"
        print "Total Mass Loss       =", 100 * (initialMass - totalMass)/initialMass, "%"
        print bold + "Gravity vs Evolution" + reset
        print "EvoError/GraError     =", evolutionError/gravityError
        print bold + "Totals" + reset
        totalGravityError = totalGravityError * (1 + gravityError/100)
        print "Total Gravity Error   =", 100 - totalGravityError * 100, "%"
        totalEvolutionError = totalEvolutionError * (1 + evolutionError/100)
        print "Total Evolution Error =", 100 - totalEvolutionError * 100, "%"
        print "Grav - Evo Error      =", (100 - totalGravityError * 100) - (100 - totalEvolutionError * 100), "%"
        
        isBound = gravity.kinetic_energy.value_in(units.J) / gravity.potential_energy.value_in(units.J)

        print "Kin / Pot             =", isBound
        
        if isBound < 0.0:
          print "Bound?                = yes"
        else:
          print "Bound?                = no"
        
        print "Evolved model to t    =", bold + str(time) + reset
        print "Evolved model to t    =", str(convert_nbody.to_nbody( time.value_in(units.Myr)| units.Myr))

        file.write(str(time.value_in(units.Myr)))
        file.write(' ')
        file.write(str(totalMass))
        file.write(' ')
        file.write(str(preEvolutionMass - totalMass))
        file.write(' ')
        file.write(str((preEvolutionMass - totalMass)/preEvolutionMass * 100))
        file.write(' ')
        file.write(str(initialMass - totalMass))
        file.write(' ')
        file.write(str(100 * (initialMass - totalMass)/initialMass))
        file.write(' ')
        file.write(str(postEvolutionEnergy))
        file.write(' ')
        file.write(str(abs(totalEnergyError)))
        file.write(' ')
        file.write(str(gravityError))
        file.write(' ')
        file.write(str(evolutionError))
        file.write(' ')
        file.write(str(100 - totalGravityError * 100))
        file.write(' ')
        file.write(str(100 - totalEvolutionError * 100))
        file.write(' ')
        file.write(str(evolutionError/gravityError))
        file.write(' ')
        """
        file.write(str(particles.kinetic_energy().value_in(units.J)))
        file.write(' ')
        file.write(str(particles.potential_energy().value_in(units.J)))
        file.write(' ')
        """
        file.write(str(isBound))
        file.write('\n')

        #Count number of binaries
        stack = []
        nBinaries = 0
        nPotentialBinaries = 0
        nHardBinaries = 0
        averageKineticEnergy = 0.0
        for stars in particles:
          stackSize = len(stack)
          distanceSquared = pow(stars.x.value_in(units.m), 2.0) + pow(stars.y.value_in(units.m), 2.0) + pow(stars.z.value_in(units.m), 2.0)
          xPosition = stars.x.value_in(units.m)
          yPosition = stars.y.value_in(units.m)
          zPosition = stars.z.value_in(units.m)
          lastMass = stars.mass.value_in(units.kg)
          xVelocity = stars.vx.value_in(units.ms)
          yVelocity = stars.vy.value_in(units.ms)
          zVelocity = stars.vz.value_in(units.ms)

          stack.append([stackSize])
          stack[stackSize].append(distanceSquared)  # x1
          stack[stackSize].append(xPosition)        # x2
          stack[stackSize].append(yPosition)        # x3
          stack[stackSize].append(zPosition)        # x4
          stack[stackSize].append(lastMass)         # x5
          stack[stackSize].append(xVelocity)        # x6
          stack[stackSize].append(yVelocity)        # x7
          stack[stackSize].append(zVelocity)        # x8

          velocitySquared = pow(stars.vx.value_in(units.ms), 2.0) + pow(stars.vy.value_in(units.ms), 2.0) + pow(stars.vz.value_in(units.ms), 2.0)
          averageKineticEnergy += 0.5 * velocitySquared * stars.mass.value_in(units.kg)

        averageKineticEnergy = averageKineticEnergy / (stackSize + 1)

        #Sort stack by distance from center
        stack.sort(key=lambda x: x[1])

        #Look at the $nCompare closest "neighbors" (because not sure if they are actually neighbors, but if they are a binary, the odds are in favor of them being neighbors in this definition)
        nCompare = 2

        for x in range (0, stackSize):
          if x + nCompare <= stackSize:
            # While there are stars to compare with, do the computations
            for nNeighbors in range (1, nCompare + 1):
              # 1. Look at the Energy of the pairs
              y = x + nNeighbors
              xSeparationSquared = pow(stack[x][2] - stack[y][2], 2.0)
              ySeparationSquared = pow(stack[x][3] - stack[y][3], 2.0)
              zSeparationSquared = pow(stack[x][4] - stack[y][4], 2.0)
              xDeltaVelocitySquared = pow(stack[x][6] - stack[y][6], 2.0)
              yDeltaVelocitySquared = pow(stack[x][7] - stack[y][7], 2.0)
              zDeltaVelocitySquared = pow(stack[x][8] - stack[y][8], 2.0)
              separation = pow(xSeparationSquared + ySeparationSquared + zSeparationSquared, 0.5)
              deltaVelocitySquared = xDeltaVelocitySquared + yDeltaVelocitySquared + zDeltaVelocitySquared
              reducedMass = stack[x][5] * stack[y][5] / (stack[x][5] + stack[y][5])
              kineticEnergy = 0.5 * reducedMass * deltaVelocitySquared
              potentialEnergy = - 6.673e-11 * stack[x][5] * stack[y][5] / separation
              totalEnergy = kineticEnergy + potentialEnergy

              if totalEnergy < 0.0:
                # Compute Semi Major Axis (a)
                # Specific orbital energy = E_tot / mu, a = - G (m1 + m2) / 2 epsilon
                orbitalEnergy = totalEnergy / reducedMass
                semiMajorAxis = - 6.673e-11 * (stack[x][5] + stack[y][5]) / (2 * orbitalEnergy)
                nPotentialBinaries += 1

                # 2. Look only at neighbors with a < 0.1 pc (Kroupa et al - "The formation of very wide binaries")
                if semiMajorAxis < 0.1 * 3.08568025e16:
                  # 3. Look at the eccentricity
                  # Calculate specific angular momentum squared (L / mu)^2 = (r x v)^2
                  angularMomentum = pow((stack[x][2] - stack[y][2]) * (stack[x][7] - stack[y][7]) - (stack[x][3] - stack[y][3]) * (stack[x][6] - stack[y][6]), 2.0)
                  angularMomentum += pow((stack[x][3] - stack[y][3]) * (stack[x][8] - stack[y][8]) - (stack[x][4] - stack[y][4]) * (stack[x][7] - stack[y][7]), 2.0)
                  angularMomentum += pow((stack[x][4] - stack[y][4]) * (stack[x][6] - stack[y][6]) - (stack[x][2] - stack[y][2]) * (stack[x][8] - stack[y][8]), 2.0)
                  eccentricity = pow(1.0 - angularMomentum / (6.673e-11 * (stack[x][5] + stack[y][5]) * semiMajorAxis), 0.5)
                  if eccentricity <= 1.0 and eccentricity >= 0.0:
                    nBinaries += 1
                    if - potentialEnergy > averageKineticEnergy:
                      nHardBinaries +=1

        nSoftBinaries = nBinaries - nHardBinaries
        print "N Possible Binaries   =", nPotentialBinaries
        print "Number of Binaries    =", nBinaries
        print "N Hard Binaries       =", nHardBinaries
        print "N Soft Binaries       =", nSoftBinaries
        binaryFile.write(str(time.value_in(units.Myr)))
        binaryFile.write(' ')
        binaryFile.write(str(nPotentialBinaries))
        binaryFile.write(' ')
        binaryFile.write(str(nBinaries))
        binaryFile.write(' ')
        binaryFile.write(str(nHardBinaries))
        binaryFile.write(' ')
        binaryFile.write(str(nSoftBinaries))
        binaryFile.write('\n')

        #Compute half-mass radius
        #Expensive calculation, do this only once every N loops (20?)
        Ncheck = 20.0
        Nloop = Nloop + 1.0
        if Nloop/Ncheck == int(Nloop/Ncheck):
          lagrangeMass = 0.0
          stack = []

          for stars in particles:
            stackSize = len(stack)

            distanceSquared = pow(stars.x.value_in(units.parsec), 2.0) + pow(stars.y.value_in(units.parsec), 2.0) + pow(stars.z.value_in(units.parsec), 2.0)
            lagrangeMass += stars.mass.value_in(units.MSun)
            lastMass = stars.mass.value_in(units.MSun)
            lastMassSI = stars.mass.value_in(units.kg)

            # Legenda:
            # stack[x][0] = ID | stack[x][1] = distance^2 | stack[x][2] = mass | stack[x][3] = Etot | stack[x][4] = x-position | stack[x][5] = y-position | stack[x][6] = z-position
            stack.append([stackSize])
            stack[stackSize].append(distanceSquared)
            stack[stackSize].append(lastMass)

            velocitySquared = pow(stars.vx.value_in(units.ms), 2.0) + pow(stars.vy.value_in(units.ms), 2.0) + pow(stars.vz.value_in(units.ms), 2.0)
            Etot = particles[stackSize].potential().value_in(units.ms**2) * stars.mass.value_in(units.kg) + 0.5 * velocitySquared * stars.mass.value_in(units.kg)
            stack[stackSize].append(Etot)
            stack[stackSize].append(stars.x.value_in(units.m))
            stack[stackSize].append(stars.y.value_in(units.m))
            stack[stackSize].append(stars.z.value_in(units.m))

          stackSize = len(stack)
          stack.sort(key=lambda x: x[1])

          # Stars sorted by distance
          # Find half mass radius
          tmpMass = 0.0
          for x in range (0, stackSize):
            tmpMass += stack[x][2]
            if tmpMass/lagrangeMass >= 0.5:
              RHalf = stack[x][1]
              print "Half-mass radius      =", RHalf, "pc"
              print "Total mass            =", lagrangeMass, "MSun"
              print "M(1/2)                =", tmpMass, "MSun"
              print "M(Tot)/M(1/2)         =", tmpMass/lagrangeMass
              break

          RHalf = pow(5.0, 2.0) * RHalf

          # Found half mass radius
          # Find stars on bigger distance than 20 * half-mass radius
          nStarsUnbound = 0.0
          nStarsOutside = 0.0
          nBound = 0.0
          realCOMx = 0.0
          realCOMy = 0.0
          realCOMz = 0.0
          realMass = 0.0
  
          for x in range (0, stackSize):
            if stack[x][1] > RHalf:
              nStarsOutside += 1.0
              if stack[x][3] > 0.0:
                nStarsUnbound += 1.0

            # Recalculate Mass
            if stack[x][1] < RHalf or stack[x][3] < 0.0:
              realMass += stars.mass.value_in(units.MSun)

          for x in range (0, stackSize):
            if stack[x][1] < RHalf or stack[x][3] < 0.0:
              #Center of Mass: R = SUM r(i) m(i) / Mtot
              realCOMx += stack[x][4] * stack[x][2] / realMass
              realCOMy += stack[x][5] * stack[x][2] / realMass
              realCOMz += stack[x][6] * stack[x][2] / realMass
              nBound += 1.0

          realCOM = [realCOMx, realCOMy, realCOMz] | units.m
          fractionUnbound = nStarsUnbound / (stackSize + 1)
          print bold + "realCOM recalculated" + reset
          unboundFile.write(str(time.value_in(units.Myr)))
          unboundFile.write(' ')
          unboundFile.write(str(nStarsUnbound))
          unboundFile.write(' ')
          unboundFile.write(str(fractionUnbound))
          unboundFile.write('\n')

        print "NStars R > 5 R_0.5    =", nStarsOutside
        print "Nr of unbound stars   =", nStarsUnbound
        print "Fraction unbound      =", fractionUnbound
        print "Nr of bound stars     =", nBound
        print "Fraction bound        =", nBound / (stackSize + 1)
        print "Real Center Of Mass   =", realCOM.value_in(units.parsec)
        print "Center of Mass        =", particles.center_of_mass().value_in(units.parsec), "pc"
        print "Center of Velocity    =", particles.center_of_mass_velocity()
     
        #Compute Lagrangian Radii
        stack = []
        for stars in particles:
          #move particles to center of mass
          stars.position = stars.position - realCOM
          distanceSquared = pow(stars.x.value_in(units.parsec), 2.0) + pow(stars.y.value_in(units.parsec), 2.0) + pow(stars.z.value_in(units.parsec), 2.0) 
          stack.append(distanceSquared)
        stackSize = len(stack)
        stack.sort()
        lagrangeFile.write(str(time.value_in(units.Myr)))
        lagrangeFile.write(' ')
        # 0.5 %
        lagrangeFile.write(str(pow(stack[int(round(len(stack)/200.0, 0))], 0.5)))
        lagrangeFile.write(' ')
        # 1 %
        lagrangeFile.write(str(pow(stack[int(round(len(stack)/100.0, 0))], 0.5)))
        lagrangeFile.write(' ')
        # 5 %
        lagrangeFile.write(str(pow(stack[int(round(len(stack)/20.0, 0))], 0.5)))
        lagrangeFile.write(' ')
        # 10 %
        lagrangeFile.write(str(pow(stack[int(round(len(stack)/10.0, 0))], 0.5)))
        lagrangeFile.write(' ')
        # 25 %
        lagrangeFile.write(str(pow(stack[int(round(len(stack)/4.0, 0))], 0.5)))
        lagrangeFile.write(' ')
        # 50 %
        lagrangeFile.write(str(pow(stack[int(round(len(stack)/2.0, 0))], 0.5)))
        lagrangeFile.write(' ')
        # 75 %
        lagrangeFile.write(str(pow(stack[int(round(len(stack)/(100.0/75.0), 0))], 0.5)))
        lagrangeFile.write(' ')
        # 90 %
        lagrangeFile.write(str(pow(stack[int(round(len(stack)/(100.0/90.0), 0))], 0.5)))
        lagrangeFile.write(' ')
        
        lagrangeFile.write('\n')
        
        # Write snapshots
        snapshotLocation = "snapshot/" + str(time.value_in(units.Myr)) + ".dat"
        regexPattern = re.compile('\.')
        snapshotLocation = regexPattern.sub('-', snapshotLocation, count=1)
        snapshotCommand = "echo \" \" > " + snapshotLocation
        plotLocation = "plot/" + str(time.value_in(units.Myr)) + ".png"
        plotLocation = regexPattern.sub('-', plotLocation, count=1)
        system(snapshotCommand)
        snapshotFile = open(snapshotLocation, 'a')
        for x in particles:
          snapshotFile.write(str(x.position.x.value_in(units.parsec)))
          snapshotFile.write(' ')
          snapshotFile.write(str(x.position.y.value_in(units.parsec)))
          snapshotFile.write(' ')
          snapshotFile.write(str(x.position.z.value_in(units.parsec)))
          snapshotFile.write('\n')          
        snapshotFile.close        
        # Plot snapshots
        comFile = open('COM.dat', 'w')
        comFile.write(str(particles.center_of_mass().x.value_in(units.parsec)) + ' ' + str(particles.center_of_mass().y.value_in(units.parsec)) + ' ' + str(particles.center_of_mass().z.value_in(units.parsec)) + '\n')
        comFile.close
        plotTitle =  str(time.value_in(units.Myr))
        plotFile = open('snap.gpl', 'w')
        plotFile.write('set title \"' + plotTitle  + ' Myr\"\nset terminal png font verdana 10 x000000 xffffff\nset output \"' + plotLocation + '\"\nunset key\nset ylabel \"y (pc)\"\nset xlabel \"x (pc)\"\nset xrange[-10:10]\nset yrange[-10:10]\nset pointsize 1\nplot \"' + snapshotLocation + '\" using 1:2 title \'\' with  points pointtype 7 linetype rgb \"yellow\", \"COM.dat\" using 1:2 title \'\' with points linewidth 6 linetype rgb \"red\"\n')
        plotFile.close
        system("gnuplot snap.gpl")
    
    if os.path.exists('small.hdf5'):
        os.remove('small.hdf5')
    storage = store.StoreHDF("small.hdf5")
    storage.store(particles)
   
    del gravity
    del stellar_evolution

    # gnuplot
    lagrangeFile.close()
    file.close()
    binaryFile.close()
    unboundFile.close()
    system("gnuplot plot.gpl > /dev/null 2&>1")
    
    #plot_particles(particles, name_of_the_figure)
    
        

        
def test_simulate_small_cluster():
    """test_simulate_small_cluster
    This method is found by the testing framework and automatically
    run with all other tests. This method simulates
    a too small cluster, this is done to limit the testing time.
    """
    assert is_mpd_running()
    simulate_small_cluster(4, 4 | units.Myr)
    
if __name__ == '__main__':
    simulate_small_cluster(int(sys.argv[1]), int(sys.argv[2]) | units.Myr)#, sys.argv[3])