from amuse.community.nbody6xx.interface import Nbody6xx
from amuse.units import nbody_system
import numpy as np
import matplotlib.pyplot as plt
from amuse.lab import *
from time import sleep

if __name__=="__main__":
    converter = nbody_system.nbody_to_si(700 | units.MSun, 1 | units.parsec)


    instance = Nbody6xx(converter,redirection="none")
    instance.initialize_code()
    instance.set_kz(20,0)
    pts = new_plummer_model(1000, converter)
    pts.radius = 1 | units.RSun
    pts.mass = new_salpeter_mass_distribution(len(pts.mass))
    pts.scale_to_standard(convert_nbody=converter)
    instance.set_rbar(1)
    instance.set_zmbar(0.7)
    instance.particles.add_particles(pts)
    #print "potential energy = ", instance.get_potential_energy()
    print "total mass = ", instance.get_total_mass()
    print "total radius = ", instance.get_total_radius()
    t_end = 0.1| units.Myr
    t=instance.get_time()
    #plt.ion()
    while (t<t_end):
        instance.evolve_model(t)
        t=instance.get_time()
        print "t = ", t.value_in(units.Myr)
        print "n = ", instance.get_number_of_particles()
        #print instance.particles.x.value_in(nbody_system.length)
        #plt.clf()
        #plt.scatter(instance.particles.x.value_in(units.AU),instance.particles.y.value_in(units.AU))
        #plt.draw()
        #plt.show()
        print instance.get_time_step()
        t += instance.get_time_step()
        sleep(0.1)

    instance.stop()
    #plt.close()
