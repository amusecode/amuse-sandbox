"""
   ./ibis-deploy.sh --gui -j lib/ibis/strw.misc.jungle 

   Nbody integration of N particles in N-body units from t=0 to
   t_end=1 Myr.  The initial cluster is a King (1966) model with
   dimension-less depth of the potential of W0=7. The initial
   distribution of stars is in virial equilibrium.  At this moment a
   4th order Hermite integrator is used for the integration.  Stellar
   masses are selected randomly from a Salpeter initial mass function
   between a minimum mass of Mmin=0.1MSun and Mmax=100MSun.  In order
   to assure proper scaling to astrophysical units, we have to define
   the cluster radius in physical units, in this case, we opted for a
   virial radius of 1pc.
   cluster in orbit in static potential
"""
import math
import sys
import numpy
from matplotlib import pyplot
from amuse.plot import loglog, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser

from amuse.couple import bridge

# import the various N-body codes
from amuse.community.hermite0.interface import Hermite
from amuse.community.fi.interface import Fi
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.bonsai.interface import Bonsai

# import stellar evolution code
from amuse.community.sse.interface import SSE

from amuse.ext.galactics_model import new_galactics_model

def main(N=10, W0=7.0, t_end=10, dt=1, filename="nbody.hdf5", Rvir=1, Mmin=0.1, Mmax=100, rgc=8500., vgc=220, z=0.02, Ngal=1000, Mgal=1.e+10, Rgal=3500):
    numpy.random.seed(1111)

    t_end = t_end | nbody_system.time
    dt = dt | t_end.unit
    Rvir = Rvir | units.parsec
    Mmin = Mmin | units.MSun
    Mmax = Mmax | units.MSun
    pos = [rgc,0,0] | units.parsec
    vel = [0,vgc,0] | units.kms
    t_orb = 2*numpy.pi*pos.length()/vel.length()
    print "t_orb=", t_orb
    Mgal = Mgal | units.MSun
    Rgal = Rgal | units.parsec

    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)

    SC_converter=nbody_system.nbody_to_si(masses.sum(),Rvir)
    GC_converter=nbody_system.nbody_to_si(Mgal,Rgal)

    bodies = new_king_model(N, W0,convert_nbody=SC_converter)
    bodies.mass = masses
    Mtot_init = masses.sum()

    # star-up stellar evolution code
    stellar = SSE(channel_type="ibis", hostname="leede")
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(bodies)
    stellar.commit_particles()

    bodies.scale_to_standard(convert_nbody=SC_converter)
    bodies.position += pos
    bodies.velocity += vel
    bodies.radius = stellar.particles.radius

    channel_from_se_to_framework = stellar.particles.new_channel_to(bodies)
    channel_from_se_to_framework.copy_attributes(["mass","radius","luminosity"])

#    cluster = Hermite(SC_converter, channel_type="ibis", hostname="local")
    cluster = Hermite(SC_converter, channel_type="ibis", hostname="koppoel")

    cluster.particles.add_particles(bodies)
#    channel_from_gravity_to_framework = cluster.particles.new_channel_to(bodies)

    galaxy_particles = new_plummer_sphere(Ngal, convert_nbody=GC_converter)
#    Nbulge = 100
#    Nhalo = 100
#    Ndisk = 100
#    N = Nbulge+Nhalo+Ndisk
#    galaxy_particles = new_galactics_model(disk_number_of_particles=Ndisk, 
#                                           bulge_number_of_particles=Nbulge, 
#                                           halo_number_of_particles=Nhalo) 
#                                           do_scale = True)
#    galaxy_particles.mass = Mgal/N
#    print "M=", galaxy_particles.mass

#    galaxy = Bonsai(GC_converter, channel_type="ibis", hostname="stoffel")
#    galaxy = Bonsai(GC_converter, channel_type="ibis", hostname="local")
#    galaxy = BHTree(GC_converter)
#    galaxy = Fi(GC_converter, channel_type="ibis", hostname="stoffel")
    galaxy = Fi(GC_converter, channel_type="ibis", hostname="thorbeckegracht")
#    galaxy = BHTree(GC_converter, channel_type="ibis", hostname="local")
    galaxy.particles.add_particles(galaxy_particles)

    write_set_to_file(bodies.savepoint(0.0|units.Myr), filename, 'hdf5')
    write_set_to_file(galaxy_particles.savepoint(0.0|units.Myr), "galaxy.hdf5", 'hdf5')
    
    gravity = bridge.Bridge()
    # bridge between star cluster and MilkyWay with one-way influence
    gravity.add_system(galaxy)
    gravity.add_system(cluster, (galaxy,) )

    channel_from_framework_to_gd = bodies.new_channel_to(gravity.particles)
    channel_from_gd_to_framework = gravity.particles.new_channel_to(bodies)
    channel_from_framework_to_gal = galaxy_particles.new_channel_to(galaxy.particles)
    channel_from_gal_to_framework = galaxy.particles.new_channel_to(galaxy_particles)

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    dt_cluster = SC_converter.to_si(1.0 | nbody_system.time)
    dt_galaxy = GC_converter.to_si(1.0/32. | nbody_system.time)
    t_end = GC_converter.to_si(t_end)
    dt = GC_converter.to_si(dt)
    print dt_cluster, dt_galaxy, dt, t_end
#    gravity.timestep = min(dt_cluster, t_orb/256.)
    gravity.timestep = dt_galaxy
    time = 0.0 | t_end.unit
    while time < t_end:
        time += dt
        gravity.evolve_model(time)
#        channel_from_gravity_to_framework.copy()

        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        stellar.evolve_model()

        channel_from_gal_to_framework.copy()
        channel_from_gd_to_framework.copy()
        channel_from_se_to_framework.copy_attributes(["mass","radius","luminosity"])
        channel_from_framework_to_gd.copy_attributes(["mass"])

        write_set_to_file(bodies.savepoint(time), filename, 'hdf5')
        write_set_to_file(galaxy_particles.savepoint(time), "galaxy.hdf5", 'hdf5')

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
    result.add_option("-r", dest="rgc", type="float",default = 8500.0,
                      help="distance to the galactic center [8500] in parsec")
    result.add_option("-v", dest="vgc", type="float",default = 220.0,
                      help="orbital velotiy around the galactic center [220] in km/s")
    result.add_option("-t", dest="t_end", type="float", default = 1.0,
                      help="end time of the simulation [1] N-body units")
    result.add_option("-W", dest="W0", type="float", default = 7.0,
                      help="Dimension-less depth of the King potential (W0) [7.0]")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    result.add_option("--Mgal", dest="Mgal", type="float",default = 1.e+10,
                      help="total mass of the Galaxy [1e+10] MSun")
    result.add_option("--Rgal", dest="Rgal", type="float",default = 3500,
                      help="size of the Galaxy [3500] pc")
    result.add_option("--Ngal", dest="Ngal", type="int",default = 1000,
                      help="Number of galaxy particles [1000]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

