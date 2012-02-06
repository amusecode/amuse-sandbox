"""
   SPH simulation of an isothermal plummer distribution of gas.
"""

import sys
import numpy
from matplotlib import pyplot
from amuse.plot import loglog, xlabel, ylabel
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser

from amuse.ext.derived_grav_systems import copycat
from amuse.ext.bridge import bridge

from amuse.ic.gasplummer import new_plummer_gas_model

from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import regular_grid_unit_cube
from amuse.ext.evrard_test import body_centered_grid_unit_cube

def main(N=1000, M=1000, R=1, t_end=10, dt=1, filename="sph.hdf5"):
    t_end = t_end | units.Myr
    dt = dt | t_end.unit
    M = M | units.MSun
    R = R | units.parsec

    converter = nbody_system.nbody_to_si(M,R)

#    gas=new_plummer_gas_model(N,convert_nbody=converter, base_grid=regular_grid_unit_cube)
#    gas.h_smooth=0. | units.parsec

    gas=molecular_cloud(targetN=N,convert_nbody=converter,
                        base_grid=body_centered_grid_unit_cube).result
    gas.h_smooth=0. | units.parsec

    mu=1.4 | units.amu
    gamma1=1.6667-1
#  print 'min Temp:', (gamma1*min(gas_parts.u)*(1.4*units.amu)/constants.kB).in_(units.K)
    print 'min Temp:', (gamma1*min(gas.u)*mu/constants.kB).in_(units.K)

    sph=Fi(converter)

    sph.parameters.use_hydro_flag=True
    sph.parameters.radiation_flag=False
    sph.parameters.self_gravity_flag=False
    sph.parameters.gamma=1
    sph.parameters.isothermal_flag=True
    sph.parameters.integrate_entropy_flag=False
    sph.parameters.timestep=dt/2.
    sph.parameters.verbosity=0 

    sph.gas_particles.add_particles(gas)

#    sph.gas_particles.u += E_inject/sph.gas_particles.total_mass()

    gravity=copycat(Fi, sph, converter)

    molecular_cloud=bridge(verbose=False)
    molecular_cloud.add_system(sph,(gravity,),False)

    time = 0.0 | t_end.unit
    while time < t_end:
        time += dt
        molecular_cloud.evolve_model(time,timestep=dt) 

        Ek=molecular_cloud.kinetic_energy.value_in(1.e51*units.erg)
        Ep=molecular_cloud.potential_energy.value_in(1.e51*units.erg)
        Eth=molecular_cloud.thermal_energy.value_in(1.e51*units.erg)
        print 't Ek Ep Eth Ef:', time,Ek,Ep,Eth,Ek+Ep+Eth

        write_set_to_file(sph.particles, filename, 'hdf5')

#    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)


    sph.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "sph.hdf5",
                      help="output filename [sph.hdf5]")
    result.add_option("-N", dest="N", type="int",default = 100,
                      help="number of stars [10]")
    result.add_option("-M", dest="M", type="float",default = 1000,
                      help="cloud mass [1000] MSun")
    result.add_option("-R", dest="R", type="float",default = 1,
                      help="cloud radius [1] pc")
    result.add_option("-t", dest="t_end", type="float", default = 1,
                      help="end time of the simulation [1] Myr")
    result.add_option("-d", dest="dt", type="float", default = 0.1,
                      help="diagnostics time step [0.1] Myr")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

