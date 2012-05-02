import numpy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

from amuse.units import generic_unit_converter
from amuse.units import units, constants
from amuse.units.generic_unit_system import *

from amuse.community.athena.interface import Athena

from cooling import evolve_internal_energy, global_mu

def set_initial_conditions(instance, converter):
        inmem = instance.grid.copy_to_memory()
        numpy.random.seed(12345)
        inmem.rho = converter.to_si(numpy.random.uniform(0.5, 1.5,
#~                loc   = 1.0, 
#~                scale = 0.5, 
                size  = instance.grid.shape
            ) | density)
        inmem.rhovx = converter.to_si(0.0 | momentum_density)
        inmem.rhovy = converter.to_si(0.0 | momentum_density)
        inmem.rhovz = converter.to_si(0.0 | momentum_density)
        inmem.energy = (5.0e11 | units.erg / units.g) * inmem.rho
        from_model_to_code = inmem.new_channel_to(instance.grid)
        from_model_to_code.copy()
    
        instance.initialize_grid()

if __name__=="__main__":
    number_of_grid_points = 40
    halfway = int((number_of_grid_points - 1) / 2)
    
    t_end = 1.0 | units.Myr
    n_steps = 100
    dt = t_end / n_steps
    
    length_unit = 10.0 | units.parsec
    mass_unit  = (1.14 | units.amu/units.cm**3) * length_unit**3 # ~ total mass
    speed_unit = 1.0 | units.kms
    
    converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(length_unit, mass_unit, speed_unit)
    hydro_code = Athena(unit_converter = converter, number_of_workers = 2)#,redirection='none')
    hydro_code.initialize_code()
    
    hydro_code.parameters.gamma = 1.001
    hydro_code.parameters.courant_number = 0.3
    hydro_code.parameters.mesh_size = (number_of_grid_points, number_of_grid_points, number_of_grid_points)
    hydro_code.parameters.length_x = 1 | length
    hydro_code.parameters.length_y = 1 | length
    hydro_code.parameters.length_z = 1 | length
    hydro_code.parameters.x_boundary_conditions = ("periodic", "periodic")
    hydro_code.parameters.y_boundary_conditions = ("periodic", "periodic")
    hydro_code.parameters.z_boundary_conditions = ("periodic", "periodic")
    hydro_code.commit_parameters()
    
    set_initial_conditions(hydro_code, converter)
    
    times = t_end * range(n_steps) / n_steps
    for i, t in enumerate(t_end * range(n_steps) / float(n_steps)):
        print t.in_(units.Myr), hydro_code.model_time.in_(units.Myr)
        print "Hydro evolve 1"
        hydro_code.evolve_model(t + dt/2.0)
        print "Cooling"
        hydro_code.grid.energy = evolve_internal_energy(
            hydro_code.grid.energy / hydro_code.grid.rho, 
            dt, 
            hydro_code.grid.rho / global_mu) * hydro_code.grid.rho
        print (global_mu / constants.kB * (hydro_code.grid.energy / hydro_code.grid.rho).amin()).in_(units.K), 
        print (global_mu / constants.kB * (hydro_code.grid.energy / hydro_code.grid.rho).amax()).in_(units.K)
        print "Hydro evolve 2"
        hydro_code.evolve_model(t + dt)
        
        pyplot.figure(figsize = (8, 8))
        hydro_code.grid[halfway].rho
        pyplot.imshow(numpy.log10(hydro_code.grid[halfway].rho.value_in(units.amu/units.cm**3)),
            vmin=0,vmax=2,origin='lower')
#~            extent=[-LL/2,LL/2,-LL/2,LL/2],vmin=0,vmax=2,origin='lower')
        pyplot.savefig('map-%6.6i.png'%i)
        pyplot.close()
    

