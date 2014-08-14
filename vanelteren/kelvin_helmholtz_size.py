"""
Runs the Kelvin-Helmholtz Instability problem in two dimensions with Athena.
"""
import numpy
import time as T
import sys

from linear_wave_amr import EvolveHydrodynamicsCodeWithAmusePeriodicBoundariesAndNCodesWithDifferentGridsSizes as AmrCode
from linear_wave_amr import Node

#from matplotlib import pyplot
from amuse.community.athena.interface import Athena
from amuse.units.generic_unit_system import *
from amuse.io import write_set_to_file
from amuse.io import ReportTable
from amuse.datamodel import Grid

NUMBER_OF_WORKERS = 4
GAMMA = 1.4
PERTUBATION_AMPLITUDE = 0.01 | speed


def new_instance_of_hydro_code(number_of_workers=1):
    result=Athena(number_of_workers = number_of_workers)
    result.parameters.gamma = GAMMA
    result.parameters.courant_number=0.8
    result.stopping_conditions.number_of_steps_detection.enable()
    return result
    


def set_parameters(instance, n, l):
    instance.parameters.mesh_size = (n,n,1)
    
    instance.parameters.length_x = l
    instance.parameters.length_y = l
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
    
def initialize_grid(grid, n):        
    vx = 0.5 | speed
    p = 2.5 | (mass / (length * time**2))
    
    halfway = n/2 - 1
    
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
        
def simulate_kelvin_helmholtz_instability(run, n, end_time, number_of_codes_per_dimension):
    
    report = ReportTable(
        "kelvin_helmholtz_run_{0}.csv".format(run), "csv", 
        attribute_types=(None, None, None,  None), 
        attribute_names=('run', 'size', 'step', 'evolve_dt', 'boundaries_dt', 'step_dt'),
        must_store_units_in_header = False
    )
    
    all_codes = []
    live_codes = []
    codes = []
    offsets = []
    is_live = []
    grid_length = 1 | length
    number_of_codes = number_of_codes_per_dimension * number_of_codes_per_dimension
    n_per_code = n / number_of_codes_per_dimension
    l_per_code = grid_length * (1.0 / number_of_codes_per_dimension)
    y = 0 *  grid_length
    dx = [1,0,0] * l_per_code
    dy = [0,1,0] * l_per_code
    offset = [0,0,0] * grid_length
    labels = []
    for j in range(number_of_codes_per_dimension):
        for i in range(number_of_codes_per_dimension):
            code = new_instance_of_hydro_code()
            codes.append(code)
            offset = (dx * i) + (dy * j)
            offsets.append(offset)
            is_live.append(True)
            live_codes.append(code)
            labels.append("{0},{1}-{2}".format(i,j,len(codes)-1))
        
    i = 0
    for j in range(number_of_codes_per_dimension):
        base = (dx * i) + (dy * j)
        offset = base - dx
        index = (number_of_codes_per_dimension - 1) + (j * number_of_codes_per_dimension)
        offsets.append(offset)
        codes.append(codes[index])
        labels.append("ghost "+labels[index])
        is_live.append(False)
        
            
    i = number_of_codes_per_dimension - 1
    for j in range(number_of_codes_per_dimension):
        base = (dx * i) + (dy * j)
        offset = base + dx
        index = (0) + (j * number_of_codes_per_dimension)
        offsets.append(offset)
        codes.append(codes[index])
        labels.append("ghost "+labels[index])
        is_live.append(False)
        
    
    j = 0
    for i in range(-1, (number_of_codes_per_dimension + 1)):
        base = (dx * i) + (dy * j)
        offset = base - dy
        index = (i % number_of_codes_per_dimension) + ((number_of_codes_per_dimension - 1)  * number_of_codes_per_dimension)
        offsets.append(offset)
        codes.append(codes[index])
        labels.append("ghost "+labels[index])
        is_live.append(False)
        
            
    j = number_of_codes_per_dimension - 1
    for i in range(-1, (number_of_codes_per_dimension + 1)):
        base = (dx * i) + (dy * j)
        offset = base + dy
        index = (i % number_of_codes_per_dimension) + (0  * number_of_codes_per_dimension)
        offsets.append(offset)
        codes.append(codes[index])
        labels.append("ghost "+labels[index])
        is_live.append(False)
            
            
    for x,y in zip(offsets,labels):
        print x,y 
    
    for x in live_codes:
        set_parameters(x, n_per_code, l_per_code)
    
    print "setup grid"
    for is_code_alive, offset, code in zip(is_live, offsets, codes):
        if not is_code_alive:
            continue
        for x in code.itergrids():
            inmem = x.copy()
            print offset
            inmem.position += offset
            clear_grid(inmem)
            initialize_grid(inmem, n)
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
    evolve_code.report = report
    evolve_code.run = run
    evolve_code.number_of_gridpoints_per_dimension = n
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
    
def save_grid(grid, index):
    write_set_to_file(grid, 'grid_{0}.hdf'.format(index))
    #pyplot.show()
    
if __name__ in ("__main__", "__plot__"):
    grids = simulate_kelvin_helmholtz_instability(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]) | time, int(sys.argv[4]))
    #for i, g in enumerate(grids):
    #    save_grid(g, i)

