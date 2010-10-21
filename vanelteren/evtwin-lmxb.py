from amuse.lab import *
import os.path

# first version for lmxb examples
# note that most parameters are given in the 'init.run_lmxb' file
# 
# this line (line 6) is important in the init.run file:
#  2.00E+00 -1.00000E+00  0.00E+00  3.50E+00  3.40E+00  3.00E-01  2.00E+00  0.00E+00    0
#   SM       DTY          AGE       PER       BMS       ECC       P1        ENC       JMX
#  SM, mass of star, not important will be overwritten by star.mass
#  DTY, don't need to set this in the file
#  AGE, same, will be 0
#  PER, Orbital period [d] 
#  BMS, total binary mass [MSun]
#  ECC, Orbital eccentricity 
#  P1,  Spin period of the primary [d]
#  ENC,  Artificial enery rate
#  JMX,  keep 0

def run():
    directory = os.path.dirname(__file__)
    print directory
    instance = EVtwin(
        init_dat_name = os.path.join(directory,'init.dat_lmxb'),
        init_run_name = os.path.join(directory,'init.run_lmxb'),
        redirection="none",
    )
    instance.initialize_code()
    instance.commit_parameters()

    particles = Particles(1)
    particles.mass = 2.0 | units.MSun
    particles.radius = 1.0 | units.RSun

    instance.particles.add_particles(particles)

    dt = 20.0 | units.Myr
    end_time = 200 | units.Myr
    t = 0.0 | units.Myr
    while t < end_time:
        t += dt
        instance.evolve_model(t)
        print instance.particles[0].mass

    instance.stop()

if __name__ == "__main__":
    run()