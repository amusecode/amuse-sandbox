from amuse.lab import *
from amuse.ext import concurrent

import sys
import time

def main(number_of_particles):
    processes = concurrent.MPIConcurrentProcesses()
    processes.init()
    
    particles = None
    if processes.is_on_root():
        particles = new_plummer_sphere(number_of_particles)
         
        ts0 = time.time()
        potential_energy0 = particles.potential_energy(G = nbody_system.G)
        ts1 = time.time()
    shared_particles = processes.share(particles)
    
    t0 = time.time()
    shared_particles.distribute()
    t1 = time.time()
    potential_energy1 = shared_particles.potential_energy(G = nbody_system.G)
    t2 = time.time()
    
    if processes.is_on_root():
        print "calculating the potential energy took:", ts1 - ts0, " seconds, SEQ"
        print "calculating the potential energy took:", t2 - t0  , " seconds, PAR"
        print "parallel is faster by factor:",  ( ts1 - ts0) /  (t2 - t0) 
        print "time for distribute:", t1 - t0,", calculate:", t2 - t1
        print "potential energy is:", potential_energy1, " SEQ"
        print "potential energy is:", potential_energy0, " PAR"
        print "error is:", (potential_energy1 - potential_energy0)/ potential_energy0
        

if __name__ == '__main__':
    main(int(sys.argv[1]))
    
