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

from amuse.ic.gasplummer import new_plummer_gas_model
from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.evrard_test import body_centered_grid_unit_cube

def main(N=10, W0=7.0, t_end=10, dt=1, filename="nbody.hdf5", Rvir=1, Mmin=0.1, Mmax=100):
    t_end = t_end | nbody_system.time
    dt = dt | t_end.unit
    Rvir = Rvir | units.parsec
    Mmin = Mmin | units.MSun
    Mmax = Mmax | units.MSun

    masses = new_salpeter_mass_distribution(N, Mmin, Mmax)
    converter=nbody_system.nbody_to_si(masses.sum(),Rvir)
    bodies = new_king_model(N, W0,convert_nbody=converter)
    bodies.mass = masses
    bodies.radius = 0 |  Rvir.unit
    bodies.scale_to_standard(convert_nbody=converter)

    gravity = Hermite(converter)
    gravity.particles.add_particles(bodies)
    
    storage = store.StoreHDF(filename,"w")
    storage.store(bodies.savepoint(t_end))

    fs_eff = 0.5
    Ngas = 100*N
    gas=new_plummer_gas_model(Ngas,convert_nbody=converter, base_grid=body_centered_grid_unit_cube)
    gas.h_smooth=0. | units.parsec
    gas.mass=gas_parts.mass*(1-sfeff)
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    bodies.move_to_center()
    gas.move_to_center()

    star_cluster_with_gas = create_system(PhiGRAPE,Gadget2,BHTree, 
                                          converter,bodies, gas)
    star_cluster_with_gas.synchronize_model()
    t_start = star_cluster_with_gas.model_time

    t_end = converter.to_si(t_end)
    dt = converter.to_si(dt)
    time = 0.0 | t_end.unit
    while time < t_end:
        sys.evolve_model(t+dt)
        sys.synchronize_model()
        time=sys.model_time
        storage.store(bodies.savepoint(time))

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        print "T=", time, "M=", bodies.mass.sum(), 
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
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
    result.add_option("-t", dest="t_end", type="float", default = 1.0,
                      help="end time of the simulation [1] N-body units")
    result.add_option("-W", dest="W0", type="float", default = 7.0,
                      help="Dimension-less depth of the King potential (W0) [7.0]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

