"""
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
from amuse.community.bhtree.interface import BHTree
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.ph4.interface import ph4

# import stellar evolution code
from amuse.community.sse.interface import SSE

class MilkyWay_galaxy(object):
    #Generic function to get gravity at a point given the potential
    def get_gravity_at_point(self, eps, x,y,z):
        phi_0 = self.get_potential_at_point(eps, x,y,z)
        grav = AdaptingVectorQuantity()
        dpos = 0.001*(x**2+y**2+z**2).sqrt()
        phi_dx = self.get_potential_at_point(0,x+dpos,y,z) - phi_0
        phi_dy = self.get_potential_at_point(0,x,y+dpos,z) - phi_0
        phi_dz = self.get_potential_at_point(0,x,y, z+dpos) - phi_0
        return phi_dx/dpos, phi_dy/dpos, phi_dz/dpos

    def disk_and_bulge_potentials(self, x,y,z, a, b, mass):
        r = (x**2+y**2).sqrt()
        return constants.G * mass /\
            (r**2 + (a + (z**2 + b**2).sqrt())**2).sqrt()

    def halo_potential(self, x,y,z, Mc=5.0E+10|units.MSun, Rc=1.0|units.kpc**2):
        r=(x**2+y**2+z**2).sqrt()
        rr = (r/Rc)
        return -constants.G * (Mc/Rc)*(0.5*numpy.log(1 +rr**2) + numpy.arctan(rr)/rr)

    #1990ApJ...348..485P
    def get_potential_at_point(self, eps, x, y, z):
        pot_disk = self.disk_and_bulge_potentials(x,y,z,
                                                  0.0|units.kpc, 
                                                  0.277|units.kpc, 
                                                  1.12E+10|units.MSun) 
        pot_bulge = self.disk_and_bulge_potentials(x,y,z, 3.7|units.kpc, 
                                                   0.20|units.kpc, 
                                                   8.07E+10|units.MSun) 
        pot_halo = self.halo_potential(x,y,z, Mc=5.0E+10|units.MSun, 
                                       Rc=6.0|units.kpc)
        return pot_disk + pot_bulge + pot_halo

def main(N=10, W0=7.0, t_end=10, dt=1, filename="nbody.hdf5", Rvir=1, Mmin=0.1, Mmax=100, rgc=8500., vgc=220, z=0.02):
    numpy.random.seed(111)
    t_end = t_end | nbody_system.time
    dt = dt | t_end.unit
    Rvir = Rvir | units.parsec
    Mmin = Mmin | units.MSun
    Mmax = Mmax | units.MSun
    pos = [rgc,0,0] | units.parsec
    vel = [0,vgc,0] | units.kms
    t_orb = 2*numpy.pi*pos.length()/vel.length()
    print "t_orb=", t_orb

    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)

    converter=nbody_system.nbody_to_si(masses.sum(),Rvir)
    bodies = new_king_model(N, W0,convert_nbody=converter)
    bodies.mass = masses
    Mtot_init = masses.sum()

    # star-up stellar evolution code
    stellar = SSE()
    stellar.parameters.metallicity = z
    stellar.particles.add_particle(bodies)
    stellar.commit_particles()

    bodies.scale_to_standard(convert_nbody=converter)
    bodies.position += pos
    bodies.velocity += vel
    bodies.radius = stellar.particles.radius

    channel_from_se_to_framework = stellar.particles.new_channel_to(bodies)
    channel_from_se_to_framework.copy_attributes(["mass","radius","luminosity"])

#    cluster_gravity = Hermite(converter)
    cluster_gravity = PhiGRAPE(converter)
#    cluster_gravity = ph4(converter)
    cluster_gravity.particles.add_particles(bodies)
    channel_from_gravity_to_framework = cluster_gravity.particles.new_channel_to(bodies)
    
    write_set_to_file(bodies.savepoint(0.0|units.Myr), filename, 'hdf5')
    
    gravity = bridge.Bridge()
    # bridge between star cluster and MilkyWay with one-way influence
    gravity.add_system(cluster_gravity, (MilkyWay_galaxy(),) )

    channel_from_framework_to_gd = bodies.new_channel_to(gravity.particles)
    channel_from_gd_to_framework = gravity.particles.new_channel_to(bodies)

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    t_end = converter.to_si(t_end)
    dt = converter.to_si(dt)
    gravity.timestep = min(dt/4., t_orb/32.)
    time = 0.0 | t_end.unit
    while time < t_end:
        time += dt
        gravity.evolve_model(time)
        channel_from_gravity_to_framework.copy()
        write_set_to_file(bodies.savepoint(time), filename, 'hdf5')

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
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

