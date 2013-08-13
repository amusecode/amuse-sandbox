from amuse.lab import *
from amuse.community.nbody6xx.interface import Nbody6xx

from amuse.units import nbody_system
from mpl_toolkits.mplot3d import Axes3D
import pylab as plt
from time import sleep

if __name__=="__main__":
    inst = Nbody6xx(redirection="none")

    inst.initialize_code()
    inst.particles.add_particles(new_king_model(1000, 3))
    inst.commit_particles()
    t=inst.get_time()
    print "total mass: ", inst.get_total_mass()

    plt.ion()
    fig = plt.figure()
    ax3D = fig.add_subplot(111, projection='3d')

    t_end = 1 |nbody_system.time
    while (t<t_end):
        print "total mass: ", inst.get_total_mass()
        inst.evolve_model(t_end)
        t=inst.get_time()
        print "t = ", t
        print "n = ", inst.get_number_of_particles()
        x_val = inst.particles.x.value_in(nbody_system.length)
        y_val = inst.particles.y.value_in(nbody_system.length)
        z_val = inst.particles.z.value_in(nbody_system.length)
        ax3D.cla()
        ax3D.scatter(x_val, y_val, z_val ,color='r', marker='o')
        #ax3D.scatter(x_val[999], y_val[999], z_val[999], s=50000, color='b', marker='o')
        plt.draw()
        plt.show()
        #sleep(1)

    inst.stop()
    plt.close()
