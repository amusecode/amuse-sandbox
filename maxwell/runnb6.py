from amuse.community.nbody6xx.interface import Nbody6xx
from amuse.units import nbody_system
import numpy as np
import matplotlib.pyplot as plt
from amuse.lab import *
from time import sleep

if __name__=="__main__":
    converter = nbody_system.nbody_to_si(1 | units.MSun, 1 | units.AU)


    instance = Nbody6xx(converter,redirection="none")
    instance.initialize_code()
    instance.particles.add_particles(new_plummer_model(1000, converter))
    #print "potential energy = ", instance.get_potential_energy()
    print "total mass = ", instance.get_total_mass()
    print "total radius = ", instance.get_total_radius()
    t_end = 0.5| units.Myr
    t=instance.get_time()
    plt.ion()
    while (t<t_end):
        instance.evolve_model(t_end)
        t=instance.get_time()
        print "t = ", t | units.Myr
        print "n = ", instance.get_number_of_particles()
        #print instance.particles.x.value_in(nbody_system.length)
        plt.clf()
        plt.scatter(instance.particles.x.value_in(units.AU),instance.particles.y.value_in(units.AU))
        plt.draw()
        plt.show()
        sleep(0.1)

    instance.stop()
    plt.close()
