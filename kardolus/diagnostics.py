# Author: Guillermo Kardolus
# Email:  kardolus@strw.leidenuniv.nl
# Usage: python diagnostics.py <nr of particles> <max time> [input_file]

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



from amuse.io import write_set_to_file
from amuse.io import read_set_from_file
from amuse.io import store
from amuse.units import nbody_system
from amuse.units import units
from amuse.community.hermite0.interface import Hermite
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPEInterface, PhiGRAPE
from amuse.community.sse.interface import SSE
from os import system

from time import time as time_code

from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution

bold = "\033[1m"
reset = "\033[0;0m"

# Create/empty datafiles
system("echo \" \" > diagnostics.dat")
system("echo \" \" > mass.dat")
system("echo \" \" > lowMass.dat")
system("echo \" \" > mediumMass.dat")
system("echo \" \" > highMass.dat")
system("echo \" \" > lagrange.dat")
system("echo \" \" > simulation.log")
system("echo \" \" > binaries.dat")
system("echo \" \" > unbound.dat")
system("echo \" \" > distribution.dat")

# Make directories
system("mkdir snapshot 2&>1")
system("mkdir binary_snapshot 2&>1")
system("mkdir plot 2&>1")
system("mkdir binary_plot 2&>1")
system("mkdir profile_dat 2&>1")
system("mkdir mass_profile_plot 2&>1")
system("mkdir density_profile_plot 2&>1")

# Remove old files
system("rm -f snapshot/* 2&>1")
system("rm -f plot/* 2&>1")
system("rm -f binary_plot/* 2&>1")
system("rm -f binary_snapshot/* 2&>1")
system("rm -f profile_dat/* 2&>1")
system("rm -f mass_profile_plot/* 2&>1")
system("rm -f density_profile_plot/* 2&>1")

# Create/empty gnuplot scripts
system("echo \" \" > snap.gpl")
system("echo \" \" > binary_snap.gpl")
system("echo \" \" > profile.gpl")

file = open('diagnostics.dat', 'a')
massFile = open('mass.dat', 'a')
lowMassFile = open('lowMass.dat', 'a')
mediumMassFile = open('mediumMass.dat', 'a')
highMassFile = open('highMass.dat', 'a')
lagrangeFile = open('lagrange.dat', 'a')
logFile = open('simulation.log', 'a')
unboundFile = open('unbound.dat', 'a')
binaryFile = open('binaries.dat', 'a')

massFile.write('lowest_mass(MSun) highest_mass(Msun) middle_mass(MSun) NStars\n')
lowMassFile.write('lowest_mass(MSun) highest_mass(Msun) middle_mass(MSun) NStars\n')
mediumMassFile.write('lowest_mass(MSun) highest_mass(Msun) middle_mass(MSun) NStars\n')
highMassFile.write('lowest_mass(MSun) highest_mass(Msun) middle_mass(MSun) NStars\n')
file.write('time(Myr) total_mass(MSun) mass_loss(MSun) mass_loss(%) total_mass_loss(MSun) total_mass_loss(%) total_energy(J) total_energy_error(%) gravity_error(%) evolution_error(%) total_gravity_error(%) total_evolution_error(%) evo_error/gra_error kinetic/potential')
file.write('\n')
lagrangeFile.write('time(Myr) 0.5(%) 1(%) 5(%) 10(%) 25(%) 50(%) 75(%) 90(%)\n')

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
    print "KE:", particles.kinetic_energy().as_quantity_in(units.J)
    print "PE:", particles.potential_energy(gravity.parameters.epsilon_squared)
    print  "center of mass:", particles.center_of_mass()
    print  "center of mass velocity:", particles.center_of_mass_velocity()

    
def simulate_small_cluster(number_of_stars, end_time = 40 | units.Myr, name_of_the_figure = "test-2.svg"):
    random.seed()
    
    salpeter_masses = new_salpeter_mass_distribution(number_of_stars)
    total_mass = salpeter_masses.sum()
    
    convert_nbody = nbody_system.nbody_to_si(total_mass, 1.0 | units.parsec)
    #convert_nbody.set_as_default()
    
    particles = new_plummer_model(number_of_stars, convert_nbody);

    gravity = PhiGRAPE(convert_nbody,mode="gpu")
    gravity.initialize_code()
    gravity.parameters.timestep_parameter=0.01
    gravity.parameters.initial_timestep_parameter=0.01
    gravity.parameters.epsilon_squared = 0.000001 | units.parsec ** 2
     
    stellar_evolution = SSE()
    stellar_evolution.initialize_module_with_default_parameters() 
    
    particles.radius = 0.0 | units.RSun
    
    #Comment next line out for equal masses
    #particles.mass = salpeter_masses
    
    if len(sys.argv) > 3:
        particles = read_set_from_file(str(sys.argv[3]), format='amuse')
        print "\n\n*** READ PARTICLES FROM", str(sys.argv[3]), "***\n\n"
    
    # Write particles to HDF5 file (.hdf) and text file
        
    if len(sys.argv) < 4:
        write_set_to_file(particles, "particles.hdf", format='amuse') 
        print "\n\n*** WROTE PARTICLES TO particles.hdf ***\n\n"
    
    gravity.particles.add_particles(particles)
    gravity.initialize_particles(0.0)
    
    from_model_to_gravity = particles.new_channel_to(gravity.particles)
    from_gravity_to_model = gravity.particles.new_channel_to(particles)
   
    time = 0.0 | units.Myr    
    particles.savepoint(time)

    move_particles_to_center_of_mass(particles)

    file.write(str(time.value_in(units.Myr)))
    file.write(' ')

    #print particles
    #write_set_to_file(particles, "particles.txt", format='txt')
    #sys.exit()

    #Count number of binaries
    stack = []
    nBinaries = 0
    nPotentialBinaries = 0
    nHardBinaries = 0
    nCloseBinaries = 0
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
        stack[stackSize].append(zVelocity)          # x8

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

                    # Look at binaries within approx 10% of gravity.parameters.epsilon_squared 
                    if pow(semiMajorAxis, 2.0) | units.m**2 < gravity.parameters.epsilon_squared * 1.1 * 1.1:
                        nCloseBinaries += 1

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
    print "N Binaries Close to e =", nCloseBinaries
    binaryFile.write('0.0')
    binaryFile.write(' ')
    binaryFile.write(str(nPotentialBinaries))
    binaryFile.write(' ')
    binaryFile.write(str(nBinaries))
    binaryFile.write(' ')
    binaryFile.write(str(nHardBinaries))
    binaryFile.write(' ')
    binaryFile.write(str(nSoftBinaries))
    binaryFile.write(' ')
    binaryFile.write(str(nCloseBinaries))
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
 
    realRHalf = pow(RHalf, 0.5)
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

    # Compute average mass
    averageMass = initialMass / nStars

    # Write to log
    logFile.write('gravity             = ')
    logFile.write(str(gravity))
    logFile.write('\n')
    logFile.write('evolution           = ')
    logFile.write(str(stellar_evolution))
    logFile.write('\n')
    logFile.write('Softening parameter = ')
    logFile.write(str(gravity.parameters.epsilon_squared.value_in(units.parsec**2)))
    logFile.write(' pc^2\n')
    logFile.write('Timestep parameter  = ')
    logFile.write(str(gravity.parameters.timestep_parameter))
    logFile.write('\n')
    logFile.write('Particle radius     = ')
    logFile.write(str(particles.radius[0]))
    logFile.write('\n')
    logFile.write('N particles         = ')
    logFile.write(str(sys.argv[1]))
    logFile.write('\n')
    logFile.write('Max time            = ')
    logFile.write(str(sys.argv[2]))
    logFile.write(' Myr\n')
    logFile.write('Average Mass        = ')
    logFile.write(str(averageMass))
    logFile.write(' MSun\n')
    logFile.write('Total Initial Mass  = ')
    logFile.write(str(initialMass))
    logFile.write(' MSun\n')
    logFile.close

    Nloop = 0.0

    # Time the while loop
    t_start = time_code()

    print "\n*** Starting the evolution of", sys.argv[1], "particles over", sys.argv[2], "Myr (max) ***\n" 

    # Use these variables to find max and min core density (used to find core collapse)
    maxRhoC = 0.0
    minRhoC = 10000.0
    isMaximum = 0
    isMinimum = 0

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

            realRHalf = pow(RHalf, 0.5)
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

        print "Half mass radius      =", realRHalf, "pc"
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
        
        #Count number of binaries and compute the central density, density, and mass profile
        #Plot the positions first
        snapshotLocation = "binary_snapshot/" + str(time.value_in(units.Myr)) + ".dat"
        regexPattern = re.compile('\.')
        snapshotLocation = regexPattern.sub('-', snapshotLocation, count=1)
        snapshotCommand = "echo \" \" > " + snapshotLocation
        plotLocation = "binary_plot/" + str(time.value_in(units.Myr)) + ".png"
        plotLocation = regexPattern.sub('-', plotLocation, count=1)
        system(snapshotCommand)
        snapshotFile = open(snapshotLocation, 'a')
        
        profileDataLocation = "profile_dat/" + str(time.value_in(units.Myr)) + ".dat"
        regexPattern = re.compile('\.')
        profileDataLocation = regexPattern.sub('-', profileDataLocation, count=1)
        profileDataCommand = "echo \" \" > " + profileDataLocation
        profilePlotLocation = "profile_plot/" + str(time.value_in(units.Myr)) + ".png"
        profilePlotLocation = regexPattern.sub('-', profilePlotLocation, count=1)
        system(profileDataCommand)
        profileDataFile = open(profileDataLocation, 'a')
        
        stack = []
        nBinaries = 0
        nPotentialBinaries = 0
        nHardBinaries = 0
        nCloseBinaries = 0
        averageKineticEnergy = 0.0
        centralMass = 0.0
        massDistribution = 0.0
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
            stack[stackSize].append(zVelocity)          # x8

            velocitySquared = pow(stars.vx.value_in(units.ms), 2.0) + pow(stars.vy.value_in(units.ms), 2.0) + pow(stars.vz.value_in(units.ms), 2.0)
            averageKineticEnergy += 0.5 * velocitySquared * stars.mass.value_in(units.kg)

        averageKineticEnergy = averageKineticEnergy / (stackSize + 1)

        #Sort stack by distance from center
        stack.sort(key=lambda x: x[1])

        #Look at the $nCompare closest "neighbors" (because not sure if they are actually neighbors, but if they are a binary, the odds are in favor of them being neighbors in this definition)
        nCompare = 2

        for x in range (0, stackSize):
            if x < 5:
                centralMass += stack[x][5]
            if x == 4:
                centralDensity = centralMass / (4.0/3.0 * 3.14159265 * pow(stack[x][1], 1.5)) * 1.4771869e19
                if Nloop == 1.0:
                    initialRhoC = centralDensity
            if Nloop/Ncheck == int(Nloop/Ncheck):
                # Compute mass and density profile once every 20 loops
                massDistribution += stack[x][5]
                specificRadius = pow(stack[x][1], 0.5)
                specificVolume = 4.0/3.0 * 3.14159265 * pow(specificRadius, 3.0)
                densityDistribution = massDistribution / specificVolume * 1.4771869e19
                profileDataFile.write(str(specificRadius/3.08568025e16))
                profileDataFile.write(' ')
                profileDataFile.write(str(massDistribution/1.98892e30))
                profileDataFile.write(' ')
                profileDataFile.write(str(densityDistribution))
                profileDataFile.write('\n')
              
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

                        # Look at binaries within gravity.parameters.epsilon_squared 
                        if pow(semiMajorAxis, 2.0) | units.m**2 < gravity.parameters.epsilon_squared * 1.1 * 1.1:
                            nCloseBinaries += 1

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
                                snapshotFile.write(str(stack[x][2]/3.08568025e16))
                                snapshotFile.write(' ')
                                snapshotFile.write(str(stack[x][3]/3.08568025e16))
                                snapshotFile.write(' ')
                                snapshotFile.write(str(stack[x][4]/3.08568025e16))
                                snapshotFile.write('\n')
                                snapshotFile.write(str(stack[y][2]/3.08568025e16))
                                snapshotFile.write(' ')
                                snapshotFile.write(str(stack[y][3]/3.08568025e16))
                                snapshotFile.write(' ')
                                snapshotFile.write(str(stack[y][4]/3.08568025e16))
                                snapshotFile.write('\n')

                                if - potentialEnergy > averageKineticEnergy:
                                    nHardBinaries += 1
                      
        snapshotFile.close
        # Plot snapshots
        plotTitle =  str(time.value_in(units.Myr))
        binaryPlotFile = open('binary_snap.gpl', 'w')
        binaryPlotFile.write('set title \"' + plotTitle  + ' Myr\"\nset terminal png font verdana 10 x000000 xffffff\nset output \"' + plotLocation + '\"\nunset key\nset ylabel \"y (pc)\"\nset xlabel \"x (pc)\"\nset autoscale fix\nset pointsize 1\nplot \"' + snapshotLocation + '\" using 1:2 title \'\' with  points pointtype 7 linetype rgb \"yellow\", \"COM.dat\" using 1:2 title \'\' with points linewidth 6 linetype rgb \"red\"\n')
        binaryPlotFile.close
        system("gnuplot binary_snap.gpl 2&>1")
        profileDataFile.close
        profilePlotTitle = str(time.value_in(units.Myr))
        profilePlotFile = open('profile.gpl', 'w')
        profilePlotFile.write('set title \"' + profilePlotTitle  + ' Myr\"\nset terminal png font verdana 10\nset output \"mass_' + profilePlotLocation + '\"\nunset key\nset ylabel \"mass (MSun)\"\nset xlabel \"radius (pc)\"\nset xrange[0:10]\nplot \"' + profileDataLocation + '\" using 1:2 title \'\' with lines\nreset\nset title \"' + profilePlotTitle  + ' Myr\"\nset terminal png font verdana 10\nset output \"density_' + profilePlotLocation + '\"\nunset key\nset ylabel \"density (MSun pc^-3)\"\nset xlabel \"radius (pc)\"\nset xrange[0:10]\nplot \"' + profileDataLocation + '\" using 1:3 title \'\' with lines\n')
        profilePlotFile.close
        system("gnuplot profile.gpl 2&>1")
        
        # Find maximum central density
        if centralDensity > 15.0 * initialRhoC and isMaximum == 0 and nBinaries > 2:
            if centralDensity > maxRhoC:
                maxRhoC = centralDensity
                collapseTime = time.value_in(units.Myr)
                collapseNbodyTime = str(convert_nbody.to_nbody( time.value_in(units.Myr)| units.Myr))
                isMinimum = 1
            
        if centralDensity < 0.5 * initialRhoC and isMinimum == 1:
            isMaximum = 1
            criticalTime = time.value_in(units.Myr)
              
        nSoftBinaries = nBinaries - nHardBinaries
        print "N Possible Binaries   =", nPotentialBinaries
        print "Number of Binaries    =", nBinaries
        print "N Hard Binaries       =", nHardBinaries
        print "N Soft Binaries       =", nSoftBinaries
        print "N Binaries close to e =", nCloseBinaries
        print "Central Density       =", centralDensity, "MSun pc^-3"
        print "Initial rho_c         =", initialRhoC, "MSun pc^-3"
        print "Central/Init Density  =", centralDensity/initialRhoC
        
        if isMinimum == 1:
            if isMaximum == 0:
                print bold + "Possible CC at t      =", collapseTime, "Myr (", collapseNbodyTime, ")"
                print "Possible CC at rho_c  =", maxRhoC, "MSun pc^-3 (", maxRhoC/initialRhoC, "rho_c(0) )" + reset
            
        if isMaximum == 1:
            print bold + "Core Collapse at t    =", collapseTime, "Myr (", collapseNbodyTime, ")"
            print "Core Collapse at rhoC =", maxRhoC, "MSun pc^-3 (", maxRhoC/initialRhoC, "rho_c(0) )"
            print "CC determined at t    =", criticalTime, "Myr" + reset
        
        binaryFile.write(str(time.value_in(units.Myr)))
        binaryFile.write(' ')
        binaryFile.write(str(nPotentialBinaries))
        binaryFile.write(' ')
        binaryFile.write(str(nBinaries))
        binaryFile.write(' ')
        binaryFile.write(str(nHardBinaries))
        binaryFile.write(' ')
        binaryFile.write(str(nSoftBinaries))
        binaryFile.write(' ')
        binaryFile.write(str(nCloseBinaries))
        binaryFile.write(' ')
        binaryFile.write(str(centralDensity))
        binaryFile.write('\n')
        
        fractionBinaries = nHardBinaries / (stackSize + 1.0)
        
        # Write snapshots
        snapshotLocation = "snapshot/" + str(time.value_in(units.Myr)) + ".dat"
        binaryLocation = "binary_snapshot/" + str(time.value_in(units.Myr)) + ".dat"
        regexPattern = re.compile('\.')
        snapshotLocation = regexPattern.sub('-', snapshotLocation, count=1)
        binaryLocation = regexPattern.sub('-', binaryLocation, count=1)
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
        plotFile.write('set title \"' + plotTitle  + ' Myr\"\nset terminal png font verdana 10 x000000 xffffff\nset output \"' + plotLocation + '\"\nunset key\nset ylabel \"y (pc)\"\nset xlabel \"x (pc)\"\nset xrange[-10:10]\nset yrange[-10:10]\nset pointsize 1\nplot \"' + snapshotLocation + '\" using 1:2 title \'\' with  points pointtype 7 linetype rgb \"yellow\", \"COM.dat\" using 1:2 title \'\' with points linewidth 6 linetype rgb \"red\", \"' + binaryLocation + '\" using 1:2 title \'\' with  points pointtype 7 linetype rgb \"blue\"\n')
        plotFile.close
        system("gnuplot snap.gpl 2&>1")

        # Exit if Core has collapsed
        if fractionBinaries > 0.01:
            time = end_time
            if isMaximum == 1:
                print "\n*** CORE COLLAPSED AT", collapseTime,"MYR AND TOO MANY BINARIES *** EXITING ***\n"
            if isMaximum == 0:
                print "\n*** TOO MANY BINARIES *** EXITING ***\n"
            print "Execution time: ", time_code() - t_start, "s (", (time_code() - t_start) / 60, "m )\n"
    
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
    
    
        

        
def test_simulate_small_cluster():
    """test_simulate_small_cluster
    This method is found by the testing framework and automatically
    run with all other tests. This method simulates
    a too small cluster, this is done to limit the testing time.
    """
    assert is_mpd_running()
    simulate_small_cluster(4, 4 | units.Myr)
    
if __name__ == '__main__':
    simulate_small_cluster(int(sys.argv[1]), int(sys.argv[2]) | units.Myr)
