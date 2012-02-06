"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
import sys
import numpy
from matplotlib import pyplot
from amuse.plot import loglog, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser

# import the various N-body codes
from amuse.community.hermite0.interface import Hermite
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPE

# import stellar evolution code
from amuse.community.sse.interface import SSE
#from amuse.community.mesa.interface import MESA
#from amuse.community.evtwin.interface import EVtwin
    
def main(N=10, W0=7.0, t_end=10, dt=1, filename="nbody.hdf5", Rvir=1, Mmin=0.1, Mmax= 100, z=0.02):
    t_end = t_end | nbody_system.time
    dt = dt | t_end.unit
    Rvir = Rvir | units.parsec
    Mmin = Mmin | units.MSun
    Mmax = Mmax | units.MSun

    # initialize cluster density model
    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)
    Mtot_init = masses.sum()
    converter=nbody_system.nbody_to_si(Mtot_init,Rvir)
    bodies = new_king_model(N, W0,convert_nbody=converter)
    bodies.mass = masses

    # star-up stellar evolution code
    stellar = SSE()
#    stellar = MESA()
#    stellar = EVtwin()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(bodies)
    stellar.commit_particles()

    channel_from_se_to_framework = stellar.particles.new_channel_to(bodies)
    channel_from_se_to_framework.copy_attributes(["mass","radius","luminosity"])

    # star up gravitational dynamics code
    bodies.radius = stellar.particles.radius
    bodies.move_to_center() # probably not needed....
    bodies.scale_to_standard(convert_nbody=converter)

    gravity = Hermite(converter)
    gravity.particles.add_particles(bodies)

    channel_from_framework_to_gd = bodies.new_channel_to(gravity.particles)
    channel_from_gd_to_framework = gravity.particles.new_channel_to(bodies)
    
    write_set_to_file(bodies, filename, 'hdf5')
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    t_end = converter.to_si(t_end)
    dt = converter.to_si(dt)
    time = 0.0 | t_end.unit
    while time < t_end:
        time += dt

        gravity.evolve_model(time)
        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        stellar.evolve_model()

        channel_from_gd_to_framework.copy()
        channel_from_se_to_framework.copy_attributes(["mass","radius","luminosity"])
        channel_from_framework_to_gd.copy_attributes(["mass"])

        write_set_to_file(bodies.savepoint(time), filename, 'hdf5')

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        dE = Etot_prev-Etot
        dE_se = Etot_prev_se-Etot
        Mtot = bodies.mass.sum()
        print "T=", time, 
        print "M=", Mtot, "(dM[SE]=", Mtot/Mtot_init, ")",
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot, 
        print "(dE[SE]=", dE_se/Etot, ")"
        Etot_init -= dE
        Etot_prev = Etot

    gravity.stop()
    stellar.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-d", dest="dt", type="float", default = 0.1,
                      help="diagnostics time step [0.1] N-body units")
    result.add_option("-f", dest="filename", default = "nbody.hdf5",
                      help="output filename [nbody.hdf5]")
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of stars [10]")
    result.add_option("-M", dest="Mmax", type="float",default = 100,
                      help="maximal stellar mass [100] MSun")
    result.add_option("-m", dest="Mmin", type="float",default = 0.1,
                      help="minimal stellar mass [0.1] MSun")
    result.add_option("-R", dest="Rvir", type="float",default = 1.0,
                      help="cluser virial radius [1] in parsec")
    result.add_option("-t", dest="t_end", type="float", default = 1.0,
                      help="end time of the simulation [1] N-body units")
    result.add_option("-W", dest="W0", type="float", default = 7.0,
                      help="Dimension-less depth of the King potential (W0) [7.0]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

