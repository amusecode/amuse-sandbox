from amuse.lab import *
from amuse.community.nbody6xx.interface import Nbody6xx
from amuse.units import nbody_system
import matplotlib.pyplot as plt
from time import sleep

if __name__=="__main__":
    inst = Nbody6xx(redirection="none")
    inst.initialize_code()
    inst.particles.add_particles(new_king_model(1000, 3))
    inst.commit_particles()
    t=inst.get_time()
    print "total mass: ", inst.get_total_mass()

    plt.ion()
    t_end = 5 |nbody_system.time
    while (t<t_end):
        print "total mass: ", inst.get_total_mass()
        inst.evolve_model(t_end)
        t=inst.get_time()
        print "t = ", t
        print "n = ", inst.get_number_of_particles()
        #print instance.particles.x.value_in(nbody_system.length)
        plt.clf()
        plt.scatter(inst.particles.x.value_in(nbody_system.length),inst.particles.y.value_in(nbody_system.length))
        plt.draw()
        plt.show()
        #sleep(0.1)

    inst.stop()
