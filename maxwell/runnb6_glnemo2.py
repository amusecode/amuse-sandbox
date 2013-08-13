from amuse.lab import *
from amuse.community.nbody6xx.interface import Nbody6xx
from amuse.units import nbody_system
from amuse.support import io
#from amuse.io.phigrape import Particles2Inp
from phigrape import Particles2Inp
from time import sleep

if __name__=="__main__":
    inst = Nbody6xx(redirection="none")

    inst.initialize_code()
    inst.particles.add_particles(new_king_model(10000, 3))
    inst.commit_particles()
    t=inst.get_time()
    print "total mass: ", inst.get_total_mass()

    #io.write_set_to_file(inst.particles, 'king_model_1000.phi','phigrape')
    t_end = 1 |nbody_system.time
    while (t<t_end):
        print "total mass: ", inst.get_total_mass()
        inst.evolve_model(t_end)
        t=inst.get_time()
        print "t = ", t
        print "n = ", inst.get_number_of_particles()
        sfname = 'snapshot/s_'+str(t.value_in(nbody_system.time))+'.phi'
        p2inp = Particles2Inp()
        p2inp.convert_to_inp(inst.particles, t.value_in(nbody_system.time),sfname)
        x_val = inst.particles.x.value_in(nbody_system.length)
        y_val = inst.particles.y.value_in(nbody_system.length)
        z_val = inst.particles.z.value_in(nbody_system.length)
        #sleep(1)

    inst.stop()
