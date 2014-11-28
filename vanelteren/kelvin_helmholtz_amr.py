"""
Runs the Kelvin-Helmholtz Instability problem in two dimensions with Athena.
"""
import numpy
import time as T
import sys

from linear_wave_amr import EvolveHydrodynamicsCodeWithAmusePeriodicBoundariesAndNCodesWithDifferentGridsSizes as AmrCode
from linear_wave_amr import Node

from matplotlib import pyplot
from amuse.community.athena.interface import Athena
from amuse.units.generic_unit_system import *
from amuse.datamodel import Grid

NUMBER_OF_WORKERS = 1
GAMMA = 1.4
DIMENSIONS_OF_MESH = (200, 200,1)
PERTUBATION_AMPLITUDE = 0.01 | speed

def new_instance_of_hydro_code(number_of_workers=1):
    result=Athena(number_of_workers = number_of_workers)
    result.parameters.gamma = GAMMA
    result.parameters.courant_number=0.8
    result.stopping_conditions.number_of_steps_detection.enable()
    return result
    


def set_parameters(instance):
    instance.parameters.mesh_size = DIMENSIONS_OF_MESH
    
    instance.parameters.length_x = 1 | length
    instance.parameters.length_y = (1.0 / NUMBER_OF_WORKERS) | length
    instance.parameters.length_z = 1 | length
    
    instance.parameters.x_boundary_conditions = ("interface","interface")
    instance.parameters.y_boundary_conditions = ("interface","interface")
    instance.parameters.z_boundary_conditions = ("periodic","periodic")
    
    
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
    
    outerregion = numpy.logical_or(grid.y <= 0.25 | length, grid.y >= 0.75 | length)
    innerregion = numpy.logical_and(grid.y > 0.25 | length, grid.y  < 0.75 | length)
    
    grid[outerregion].rho = 1  | density
    grid[outerregion].rhovx =  vx * grid[outerregion].rho
    
    grid[innerregion].rho = 2.0  | density
    grid[innerregion].rhovx = -vx * grid[innerregion].rho
    
    grid.energy = p / (GAMMA - 1)
        
def pertubate_grid(grid):
    grid.rhovx += grid.rho * PERTUBATION_AMPLITUDE * (numpy.random.rand(*grid.shape) - 0.5)
    grid.rhovy += grid.rho * PERTUBATION_AMPLITUDE * (numpy.random.rand(*grid.shape) - 0.5)
    
    grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
def simulate_kelvin_helmholtz_instability(end_time):
    number_of_codes = NUMBER_OF_WORKERS
    all_codes = []
    live_codes = []
    codes = []
    offsets = []
    is_live = []
    grid_length = 1 | length
    
    y = 0 *  grid_length
    dx = [1,0,0] * (grid_length)
    dy = [0,1,0] * (grid_length * 1.0 / number_of_codes)
    offset = [0,0,0] * grid_length
    
    for i in range(number_of_codes):
        code = new_instance_of_hydro_code()
        codes.append(code)
        offsets.append(offset)
        is_live.append(True)
        offset += dy
        
        live_codes.append(code)
        
    
    #periodical in y, before 
    offsets.append(offsets[0] - dy)
    codes.append(codes[number_of_codes-1])
    is_live.append(False)
    
    #periodical in y, after 
    offsets.append(offsets[number_of_codes-1] + dy)
    codes.append(codes[0])
    is_live.append(False)
    
    #periodical in x, one code length
    for i in range(number_of_codes+2):
        code = codes[i]
        offset = offsets[i]
        offsets.append(offset - dx)
        codes.append(code)
        is_live.append(False)
        
    #periodical in x, one code length
    for i in range(number_of_codes+2):
        code = codes[i]
        offset = offsets[i]
        offsets.append(offset + dx)
        codes.append(code)
        is_live.append(False)
        
    
    for x in live_codes:
        set_parameters(x)
    
    print "setup grid"
    for is_code_alive, offset, code in zip(is_live, offsets, codes):
        if not is_code_alive:
            continue
        for x in code.itergrids():
            inmem = x.copy()
            print offset
            inmem.position += offset
            clear_grid(inmem)
            initialize_grid(inmem)
            pertubate_grid(inmem)
            
            from_model_to_code = inmem.new_channel_to(x)
            from_model_to_code.copy()
            
            
    nodes = []
    for is_code_alive, offset, code in zip(is_live, offsets, codes):
        node = Node(code, offset, not is_code_alive)
        nodes.append(node)
        
    
    for x in live_codes:
        x.set_timestep(0.008 | end_time.unit)
    evolve_code = AmrCode(nodes)
    
    for x in live_codes:
        x.initialize_grid()

    t0 = T.time()
    evolve_code.init_channels()
    print "start evolve"
    if 1:
        dt = end_time / 1.0
        t = dt
        while t <= end_time:
            evolve_code.evolve_model(t, end_time)
            
            print "time : ", t
            t += dt
    dt = T.time() - t0
    print "evolve took: ", dt
        
    print "copying results"
    result = []
    
    for is_code_alive, offset, code in zip(is_live, offsets, codes):
        if not is_code_alive:
            continue
        for x in code.itergrids():
            inmem = x.copy()
            inmem.position += offset
            result.append(inmem)
            
    print "terminating code"
    
    for x in live_codes:
        x.stop()

    return result
    
def plot_grid(grid, index):
    rho = grid.rho[...,...,0].value_in(density)
    figure = pyplot.figure(figsize=(6,6))
    plot = figure.add_subplot(1,1,1)
    plot.imshow(rho, origin = 'lower')
    figure.savefig('kelvin_helmholtz_{0}.png'.format(index))
    #pyplot.show()
    
if __name__ in ("__main__", "__plot__"):
    grids = simulate_kelvin_helmholtz_instability(float(sys.argv[1]) | time)
    for i, g in enumerate(grids):
        plot_grid(g, i)
