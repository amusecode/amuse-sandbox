"""
Runs the Kelvin-Helmholtz Instability problem in two dimensions with Athena.
"""
import numpy
import sys
import time as T
from amuse.rfi.channel import AsyncRequestsPool
from amuse.community.athena.interface import Athena
from amuse.units.generic_unit_system import *
from amuse.datamodel import Grid

GAMMA = 1.4
DIMENSIONS_OF_MESH = (400,400,1)
PERTUBATION_AMPLITUDE = 0.01 | speed

def new_instance_of_hydro_code(number_of_workers=4):
    result=Athena(number_of_workers = number_of_workers)
    result.parameters.gamma = GAMMA
    result.parameters.courant_number=0.8
    result.stopping_conditions.number_of_steps_detection.enable()
    return result

def set_parameters(instance):
    instance.parameters.mesh_size = DIMENSIONS_OF_MESH
    
    instance.parameters.length_x = 1 | length
    instance.parameters.length_y = 1 | length
    instance.parameters.length_z = 1 | length
    
    instance.parameters.x_boundary_conditions = ("periodic","periodic")
    instance.parameters.y_boundary_conditions = ("periodic","periodic")
    instance.parameters.z_boundary_conditions = ("periodic","periodic")
    
    
def new_grid():
    grid = Grid.create(DIMENSIONS_OF_MESH, [1,1,1] | length)
    self.clear_grid(grid)
    return grid
    
def clear_grid(grid):
    density = mass / length**3
    momentum =  speed * density
    energy =  mass / (time**2 * length)

    grid.rho =  0.0 | density
    grid.rhovx = 0.0 | momentum
    grid.rhovy = 0.0 | momentum
    grid.rhovz = 0.0 | momentum
    grid.energy = 0.0 | energy

    return grid
    
def initialize_grid(grid):        
    vx = 0.5 | speed
    p = 2.5 | (mass / (length * time**2))
    
    halfway = DIMENSIONS_OF_MESH[0]/2 - 1
    
    outerregion = numpy.logical_or(grid.y <= 0.25 | length, grid.y >= 0.75 | length)
    innerregion = numpy.logical_and(grid.y > 0.25 | length, grid.y < 0.75 | length)
    
    grid[outerregion].rho = 1  | density
    grid[outerregion].rhovx =  vx * grid[outerregion].rho
    
    grid[innerregion].rho = 2.0  | density
    grid[innerregion].rhovx = -vx * grid[innerregion].rho
    
    grid.energy = p / (GAMMA - 1)
        
def pertubate_grid(grid):
    grid.rhovx += grid.rho * PERTUBATION_AMPLITUDE * (numpy.random.rand(*grid.shape) - 0.5)
    grid.rhovy += grid.rho * PERTUBATION_AMPLITUDE * (numpy.random.rand(*grid.shape) - 0.5)
    
    grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho

def evolve_model(code, time, endtime):
    t_unit = code.get_timestep().unit
    while ((code.model_time - time)/time) < -1e-14:
        timesteps = [] | t_unit
        timesteps.append(code.get_timestep())
        min_timestep = timesteps.min()
        print "timestep:", min_timestep
        #if code.model_time == 0.0 | t_unit:
        #if code.model_time == 0.0 | t_unit:
        #    min_timestep = time
            
        if (min_timestep + code.model_time) >= time and time == endtime:
            code.parameters.must_evolve_to_exact_time = True
        #print min_timestep
        code.set_timestep(min_timestep)
        
        t0 = T.time()
        pool = AsyncRequestsPool()
        request = code.evolve_model.async(time)
        pool.add_request(request, lambda x : x.result())
        pool.waitall()
        dt = T.time() - t0
        print "evolve model:", dt
            
        print "===="
        print code.model_time, (code.model_time - time)/time
        print "===="
        
def simulate_kelvin_helmholtz_instability(steps_per_evolve, total_steps):
    instance=new_instance_of_hydro_code()
    set_parameters(instance)
    instance.parameters.stopping_conditions_number_of_steps = steps_per_evolve
    
    print "setup grid"
    for x in instance.itergrids():
        inmem = x.copy()
        
        clear_grid(inmem)
        initialize_grid(inmem)
        pertubate_grid(inmem)
        
        from_model_to_code = inmem.new_channel_to(x)
        from_model_to_code.copy()
    
    print "start evolve"
    
    t0 = T.time()
    t = 10 | time
    n = 0
    
    while n < total_steps:
        t00 = T.time()
        instance.evolve_model(t)
        n += steps_per_evolve
        dt = T.time() - t00
        print "substep time:", instance.model_time
        print "evolve substep", dt
    
    dt = T.time() - t0
    print "evolve total", dt
    
    print "copying results"
    result = []
    for x in instance.itergrids():
        result.append(x.copy())

    print "terminating code"
    instance.stop()

    return result
    
    
if __name__ in ("__main__", "__plot__"):
    simulate_kelvin_helmholtz_instability(
        int(sys.argv[1]),
        int(sys.argv[2])
    )
