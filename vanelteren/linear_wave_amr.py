"""
In this script we simulate a 1d linear wave in a 2d field, 
the periodic boundaries are not handled in the code but by amuse
(much slower at time of writing)
"""
import numpy
import math
import time as T
from itertools import groupby

from amuse.test.amusetest import TestCase
from amuse.rfi.channel import AsyncRequestsPool
from amuse.rfi.channel import ASyncRequestSequence

from amuse.support.core import late
from amuse.units.quantities import VectorQuantity
from amuse.units.quantities import concatenate
from amuse.units.generic_unit_system import *
from amuse.units import generic_unit_system
from amuse.units import quantities
from amuse.units.quantities import is_quantity

from amuse.community.capreole.interface import Capreole
from amuse import io
from amuse.io import text

from amuse.datamodel import Grid
from amuse.datamodel import indexing
from amuse.datamodel import grid_attributes
from amuse.datamodel.grids import SamplePointsOnMultipleGrids
from amuse.datamodel.grids import SamplePointWithIntepolation
from amuse.datamodel.grids import SamplePointOnCellCenter
try:
    from amuse import plot
    from matplotlib import pyplot
    from matplotlib import animation
    IS_PLOT_AVAILABLE = True
except ImportError:
    IS_PLOT_AVAILABLE = False

USE_BOUNDARIES = True
        
def interpolate(value, cx, cy, cz, dx, dy, dz, ix, iy, iz, nx, ny, nz):
    iz1 = numpy.minimum(iz+1, nz-1)
    iy1 = numpy.minimum(iy+1, ny-1)
    ix1 = numpy.minimum(ix+1, nx-1)
    value000 = value[ix , iy , iz ]
    value001 = value[ix , iy , iz1]
    value010 = value[ix , iy1, iz ]
    value011 = value[ix , iy1, iz1]
    value100 = value[ix1, iy , iz ]
    value101 = value[ix1, iy , iz1]
    value110 = value[ix1, iy1, iz ]
    value111 = value[ix1, iy1, iz1]
    
    dx1 = dx / cx
    dx0 = 1 - dx1
    if nz == 1:
        dy1 = 0
        dy0 = 1
    else:
        dy1 = dy / cy
        dy0 = 1 - dy1
    if nz == 1:
        dz1 = 0
        dz0 = 1
    else:
        dz1 = dz / cz
        dz0 = 1 - dz1
    value000 *= (dx0 * dy0 * dz0)
    value001 *= (dx0 * dy0 * dz1)
    value010 *= (dx0 * dy1 * dz0)
    value011 *= (dx0 * dy1 * dz1)
    value100 *= (dx1 * dy0 * dz0)
    value101 *= (dx1 * dy0 * dz1)
    value110 *= (dx1 * dy1 * dz0)
    value111 *= (dx1 * dy1 * dz1)
    return (
        value000 + value001 + value010 + value011 +
        value100 + value101 + value110 + value111
    )
class NoCode(object):
    def __init__(self, grid):
        self.grid = grid
        
    def get_hydro_state_at_point(self, x, y, z):
        cx, cy, cz  = self.grid.cellsize()
        min  = self.grid.get_minimum_position()
        shape = self.grid.shape
        n = numpy.array([1,1,1])
        n = numpy.arange(3)[(numpy.asarray(shape) == n)]
        delta_position_x = x - self.grid[0,0,0].x
        delta_position_y = y - self.grid[0,0,0].y
        delta_position_z = z - self.grid[0,0,0].z
        index_x = numpy.asarray(numpy.floor(delta_position_x / cx), dtype=numpy.int32)
        if shape[2] == 1:
            index_z = 0 * index_x
        else:
            index_z = numpy.asarray(numpy.floor(delta_position_z / cz), dtype=numpy.int32)
        if shape[1] == 1:
            index_y = 0 * index_x
        else:
            index_y = numpy.asarray(numpy.floor(delta_position_y / cy), dtype=numpy.int32)
        dx = x - self.grid.x[index_x,index_y,index_z]
        dy = y - self.grid.y[index_x,index_y,index_z]
        dz = z - self.grid.z[index_x,index_y,index_z]
        rho = interpolate(self.grid.rho, cx, cy, cz, dx, dy, dz, index_x, index_y, index_z, *shape)
        rhovx = interpolate(self.grid.rhovx, cx, cy, cz, dx, dy, dz, index_x, index_y, index_z, *shape)
        rhovy = interpolate(self.grid.rhovy, cx, cy, cz, dx, dy, dz, index_x, index_y, index_z, *shape)
        rhovx = interpolate(self.grid.rhovz, cx, cy, cz, dx, dy, dz, index_x, index_y, index_z, *shape)
        energy = interpolate(self.grid.energy, cx, cy, cz, dx, dy, dz, index_x, index_y, index_z, *shape)
        return rho, rhovx, rhovy, rhovx, energy
    
class HydroStateAtPointTests(TestCase):
    def test1(self):
        import warnings
        warnings.simplefilter('error')
        input = Grid.create((10,1,1), [20,0,0] | generic_unit_system.length)
        input.rho = input.x * 0.1
        input.rhovx = input.x * 0.2
        input.rhovy = input.x * 0.3
        input.rhovz = input.x * 0.4
        input.energy = input.x * 0.5
        print input.x
        rho, rhovx, rhovy, rhovz, energy = NoCode(input).get_hydro_state_at_point(*([2,0,0] | generic_unit_system.length))
        self.assertAlmostRelativeEquals(rho, 0.1 * 2  | generic_unit_system.length)
        self.assertAlmostRelativeEquals(rhovz, 0.4 * 2  | generic_unit_system.length)
        
    def test2(self):
        input = Grid.create((10,1,1), [20,0,0] | generic_unit_system.length)
        input.rho = input.x * 0.1
        input.rhovx = input.x * 0.2
        input.rhovy = input.x * 0.3
        input.rhovz = input.x * 0.4
        input.energy = input.x * 0.5
        print input.x
        rho, rhovx, rhovy, rhovz, energy = NoCode(input).get_hydro_state_at_point(
            [2,4,5] | generic_unit_system.length,
            [0,0,0] | generic_unit_system.length,
            [0,0,0] | generic_unit_system.length,
        )
        self.assertAlmostRelativeEquals(rho, [0.1 * 2, 0.1 * 4, 0.1 * 5]  | generic_unit_system.length)
        
    def test3(self):
        input = Grid.create((10,20,1), [20,10,0] | generic_unit_system.length)
        input.rho = input.y * 0.1
        input.rhovx = input.x * 0.2
        input.rhovy = input.x * 0.3
        input.rhovz = input.x * 0.4
        input.energy = input.x * 0.5
        print input.x
        rho, rhovx, rhovy, rhovz, energy = NoCode(input).get_hydro_state_at_point(
            input.x,
            input.y,
            input.z
        )
        self.assertAlmostRelativeEquals(rho, input.rho)
        
def refine_grid(code, grid, offset, factor = 2):
    nx, ny, nz = grid.shape
    maxx, maxy, maxz = grid.get_maximum_position()
    minx, miny, minz = grid.get_minimum_position()
    
    nnx = min(factor, nx)
    dx = (maxx - minx) / nnx
    nny = min(factor, ny)
    dy = (maxy - miny) / nny
    nnz = min(factor, nz)
    dz = (maxz - minz) / nnz
    grids = []
    for i in range(nnx):
        x = i * dx
        for j in range(nny):
            y = i * dy
            for k in range(nnz):
                z = i *dz
                refined_grid = Grid.create(grid.shape, (dx,dy,dz,))
                rho, rhovx, rhovy, rhovz, energy = code.get_hydro_state_at_point(
                    (refined_grid.x + x).flatten(),
                    (refined_grid.y + y).flatten(),
                    (refined_grid.z + z).flatten()
                )
                refined_grid.rho = rho.reshape(refined_grid.shape)
                refined_grid.rhovx = rhovx.reshape(refined_grid.shape)
                refined_grid.rhovy = rhovy.reshape(refined_grid.shape)
                refined_grid.rhovz = rhovz.reshape(refined_grid.shape)
                refined_grid.energy = energy.reshape(refined_grid.shape)
                grids.append([quantities.as_vector_quantity([x,y,z]), refined_grid])
    return grids

class RefineGridTests(TestCase):
    def test1(self):
        input = Grid.create((10,1,1), [6,6,6] | generic_unit_system.length)
        input.rho = input.x * 0.1
        input.rhovx = input.x * 0.2
        input.rhovy = input.x * 0.3
        input.rhovz = input.x * 0.4
        input.energy = input.x * 0.5
        output = refine_grid(NoCode(input), input, [0,0,0] | generic_unit_system.length)
        self.assertEquals(len(output), 2)
        
    def test2(self):
        input = Grid.create((10,10,1), [6,6,6] | generic_unit_system.length)
        input.rho = input.y * 0.1
        input.rhovx = input.x * 0.2
        input.rhovy = input.x * 0.3
        input.rhovz = input.x * 0.4
        input.energy = input.x * 0.5
        output = refine_grid(NoCode(input), input, [0,0,0] | generic_unit_system.length)
        self.assertEquals(len(output), 4)
        
class HydroStateChannel(object):
    
    def __init__(self, from_code, to_code, to_grid,  offset, factor = 1):
        self.from_code = from_code
        self.to_code = to_code
        self.to_grid = to_grid
        self.offset = offset
        self.factor = factor
        self.prepare()
        
    def prepare(self):
        
        positions = self.to_grid.position + self.offset
        x = positions.x
        y = positions.y
        z = positions.z
        
        self.shape = x.shape
        
        self.x = x.flatten()
        self.y = y.flatten()
        self.z = z.flatten()
        self.copy_grid = self.to_grid.copy()
        
        self.to_channel = self.copy_grid.new_channel_to(self.to_grid)
        
        
    def copy(self):
        shape = self.shape
        rho, rhovx, rhovy, rhovz, rhoe = self.from_code.get_hydro_state_at_point(self.x, self.y, self.z)
       
            
        rho = rho.reshape(shape) 
        rhovx = rhovx.reshape(shape)
        rhovy = rhovy.reshape(shape)
        rhovz = rhovz.reshape(shape)
        rhoe = rhoe.reshape(shape)
        
        self.copy_grid.rho = rho
        self.copy_grid.rhovx = rhovx
        self.copy_grid.rhovy = rhovy
        self.copy_grid.rhovz = rhovz
        self.copy_grid.energy = rhoe
        
        self.to_channel.copy()
        
    def get(self, pool = None):
        request = self.from_code.get_hydro_state_at_point.async(self.x, self.y, self.z)
        pool.add_request(request, self.handle_get_result)
        
    def set(self):
        self.to_channel.copy()
        
    def handle_get_result(self, request):
        
        shape = self.shape
        rho, rhovx, rhovy, rhovz, rhoe = request.result()
        rho = rho.reshape(shape) 
        rhovx = rhovx.reshape(shape)
        rhovy = rhovy.reshape(shape)
        rhovz = rhovz.reshape(shape)
        rhoe = rhoe.reshape(shape)
        
        self.copy_grid.rho = rho
        self.copy_grid.rhovx = rhovx
        self.copy_grid.rhovy = rhovy
        self.copy_grid.rhovz = rhovz
        self.copy_grid.energy = rhoe
        
    def handle_get_result2(self, result, offset):
        
        shape = self.shape
        rho, rhovx, rhovy, rhovz, rhoe = result
        
        rho = rho[offset:offset+len(self.x)].reshape(shape) 
        rhovx = rhovx[offset:offset+len(self.x)].reshape(shape)
        rhovy = rhovy[offset:offset+len(self.x)].reshape(shape)
        rhovz = rhovz[offset:offset+len(self.x)].reshape(shape)
        rhoe = rhoe[offset:offset+len(self.x)].reshape(shape)
        
        self.copy_grid.rho = rho
        self.copy_grid.rhovx = rhovx
        self.copy_grid.rhovy = rhovy
        self.copy_grid.rhovz = rhovz
        self.copy_grid.energy = rhoe

class Node(object):
    def __init__(self, code, offset, is_ghost = False):
        self.offset = offset
        self.code = code
        self.is_ghost = is_ghost
        
    def is_live(self):
        return not self.is_ghost
    
    def setup(self):
        self.size = quantities.as_vector_quantity([
            self.code.parameters.length_x,
            self.code.parameters.length_y,
            self.code.parameters.length_z
        ])
        self.number_of_grid_points = numpy.asarray([
            self.code.parameters.nx,
            self.code.parameters.ny,
            self.code.parameters.nz
        ])
        self.cellsize = self.size / self.number_of_grid_points
        
    def contains(self, points):
        return numpy.logical_and(
            numpy.all(points >= self.offset, axis=len(points.shape)-1),
            numpy.all(points < (self.size + self.offset), axis=len(points.shape)-1)
        )
        
class EvolveHydrodynamicsCodeWithAmusePeriodicBoundariesAndNCodesWithDifferentGridsSizes(object):
    
    """
        first version, grids connect on the y axis, whole system is periodic
    """
    def __init__(self, nodes):
        self.nodes = nodes
        self.min_timestep = None
        self.channels = []
    
    def set_parameters(self):
        for node in self.nodes:
            if node.is_ghost:
                continue
            code = node.code
            code.parameters.x_boundary_conditions = ("interface","interface")
            code.parameters.y_boundary_conditions = ("interface","interface")
            code.stopping_conditions.number_of_steps_detection.enable()
            
    def get_boundaries(self, code):
        x1, x2 = code.parameters.x_boundary_conditions
        y1, y2 = code.parameters.y_boundary_conditions
        z1, z2 = code.parameters.z_boundary_conditions
        result = []
        if x1 == "interface":
            result.append(code.get_boundary_grid('xbound1'))
        if x2 == "interface":
            result.append(code.get_boundary_grid('xbound2'))
        if y1 == "interface":
            result.append(code.get_boundary_grid('ybound1'))
        if y2 == "interface":
            result.append(code.get_boundary_grid('ybound2'))
        if z1 == "interface":
            result.append(code.get_boundary_grid('zbound1'))
        if z2 == "interface":
            result.append(code.get_boundary_grid('zbound2'))
        return result
        
    def init_channels(self):
        for x in self.nodes:
            x.setup()
            
        channels = []
        for index in range(0, len(self.nodes)):
            node = self.nodes[index]
            instance = node.code
            offset = node.offset
            if node.is_ghost:
                continue
            i = 0
            for boundary in self.get_boundaries(instance):
                position = boundary.position + offset
                for other_node in self.nodes:
                    other_code = other_node.code
                    other_offset = other_node.offset
                    relative_position = position
                    selection = other_node.contains(position) 
                    i+=1
                    if not numpy.any(selection):
                        continue
                    factor = numpy.prod(node.cellsize / other_node.cellsize)
                    #print factor
                    #print offset - other_offset
                    #print i, id(instance), id(other_code)
                    channel = HydroStateChannel(other_code, node.code, boundary[selection], offset - other_offset, factor)
                    channels.append(channel)
            
            print index, len(channels)
                        
        self.channels = list(reversed(sorted(channels, key = lambda x : x.factor)))
        for x in self.channels:
            print x.factor
            
        
        self.channels_grouped_by_factor_and_code = []
        self.channels_grouped_by_factor_and_to_code = []
        for factor, group in groupby(self.channels, lambda x : x.factor): 
            grouped_by_factor = list(group)
            print "factor:", factor, len(grouped_by_factor)
            grouped_by_factor = list(sorted(grouped_by_factor, key = lambda x : id(x.from_code)))
            for from_code_id, code_group in groupby(grouped_by_factor, lambda x : id(x.from_code)): 
                grouped_by_factor_and_code = list(code_group)
                print "from code:", factor, from_code_id, len(grouped_by_factor_and_code)
                self.channels_grouped_by_factor_and_code.append(grouped_by_factor_and_code)
                
                
            grouped_by_factor = list(sorted(grouped_by_factor, key = lambda x : id(x.to_code)))
            for to_code_id, code_group in groupby(grouped_by_factor, lambda x : id(x.to_code)): 
                grouped_by_factor_and_to_code = list(code_group)
                print "to code:", factor, to_code_id, len(grouped_by_factor_and_to_code)
                grouped_by_factor_and_to_code = list(sorted(grouped_by_factor_and_to_code, key = lambda x : x.copy_grid.size))
                self.channels_grouped_by_factor_and_to_code.append(grouped_by_factor_and_to_code)
        t0 = T.time()
        for i in range(10):
            self.copy_to_boundary_cells()
        dt = T.time() - t0
        print "copy to boundary dt:", dt
        
    def copy_to_boundary_cells(self):
        if 1:
            pool = AsyncRequestsPool()
            def handle_result(request, group, offsets):
                result = request.result()
                for y, o in zip(group, offsets):
                    y.handle_get_result2(result, o)
                    
            for x in self.channels_grouped_by_factor_and_code:
                allx = []
                allz = []
                ally = []
                offsets = []
                offset = 0
                for y in x:
                    offsets.append(offset)
                    offset += len(y.x)
                    allx.append(y.x)
                    ally.append(y.y)
                    allz.append(y.z)
                allx = concatenate(allx)
                ally = concatenate(ally)
                allz = concatenate(allz)
                
                request = x[0].from_code.get_hydro_state_at_point.async(allx, ally, allz)
                pool.add_request(request, handle_result, [x, list(offsets)])
            if len(pool) > 0:
                pool.waitall()
            
            index = 0
            if 0:
                while True:
                    pool = AsyncRequestsPool()
                    for x in self.channels_grouped_by_factor_and_to_code:
                        if index < len(x):
                            y = x[index]
                            names_to_copy = list(y.to_channel.get_overlapping_attributes())
                            values = y.to_channel.get_values(names_to_copy)
                            request = y.to_channel.target.set_values_in_store_async(y.to_channel.index, names_to_copy, values)
                            pool.add_request(request, lambda z : z.result())
                    index += 1
                    if len(pool) == 0:
                        break
                    pool.waitall()
            else:
                pool = AsyncRequestsPool()
                for x in self.channels_grouped_by_factor_and_to_code:
                    values = []
                    indices = []
                    atss = []
                    names_to_copy = list(x[0].to_channel.get_overlapping_attributes())
                    for y in x:
                        ats = y.to_channel.target._private.grid._private.attribute_storage
                        atss.append(ats)
                        values.append(y.to_channel.get_values(names_to_copy))
                        i = indexing.combine_indices(y.to_channel.target._private.indices, y.to_channel.index)
                        indices.append(i)
                     
                    request = set_values_in_store_multi_async(atss, indices, names_to_copy, values)
                    pool.add_request(request, lambda z : z.result())
                if len(pool) > 0:
                    pool.waitall()
        else:
            for channel in self.channels:
                channel.copy()
        
    def evolve_model(self, time, endtime = None):
        code = self.nodes[0].code
        t_unit = code.get_timestep().unit
        step = 0
        while ((code.model_time - time)/time) < -1e-14:
            t00 = T.time()
            timesteps = [] | t_unit
            for x in self.nodes:
                if not x.is_ghost:
                    timesteps.append(x.code.get_timestep())
                    #print x.code.get_timestep()
            #print timesteps
            min_timestep = timesteps.min()
            print min_timestep, (min_timestep + code.model_time) >= time, code.model_time - time
            print (code.model_time - time)/time
            #if code.model_time == 0.0 | t_unit:
            #if code.model_time == 0.0 | t_unit:
            #    min_timestep = time
                
            if (min_timestep + code.model_time) >= time and time == endtime:
                for x in self.nodes:
                    if not x.is_ghost:
                        x.code.parameters.must_evolve_to_exact_time = True
            #print min_timestep
            for x in self.nodes:
                if not x.is_ghost:
                    x.code.set_timestep(min_timestep)
            
            t0 = T.time()
            pool = AsyncRequestsPool()
            for x in self.nodes:
                if not x.is_ghost:
                    request = x.code.evolve_model.async(time)
                    pool.add_request(request, lambda x : x.result())
            pool.waitall()
            dt = T.time() - t0
            print "evolve model:", dt
            evolve_dt = dt
                
            print "===="
            for x in self.nodes:
                if not x.is_ghost:
                    print x.code.model_time, (x.code.model_time - time)/time
            print "===="
            t0 = T.time()
            self.copy_to_boundary_cells()
            dt = T.time() - t0
            boundaries_dt = dt
            print "copy boundary values:", dt
            dtt = T.time() - t00
            print "step:", dtt
            step_dt = dtt
            step += 1
            #self.report.add_row(self.run, self.number_of_gridpoints_per_dimension, step, evolve_dt, boundaries_dt, step_dt)



def set_values_in_store_multi_async(atss, list_of_indices, attributes, list_of_quantities):
        one_dimensional_array_of_indices = None
        original_attributes = attributes
        attributes = list(attributes)
        extra_attributes = list(atss[0].extra_keyword_arguments_for_getters_and_setters)
        extra_kw = {}
        for x in extra_attributes:
            extra_kw[x] = []
        one_dimensional_values = [list() for x in attributes]
        for ats, indices, quantities in zip(atss, list_of_indices, list_of_quantities):
            array_of_indices = ats._to_arrays_of_indices(indices)
            if one_dimensional_array_of_indices is None:
                one_dimensional_array_of_indices = [list() for x in array_of_indices]
    
            one_d_v = [(x.reshape(-1) if is_quantity(x) else numpy.asanyarray(x).reshape(-1)) for x in quantities]
            extra = ats.extra_keyword_arguments_for_getters_and_setters
            for key in extra_attributes:
                value = extra[key]
                value = numpy.ones(len(one_d_v[0])) * value
                extra_kw[key].append(value)
                
            for l, x in zip(one_dimensional_values, one_d_v):
                l.append(x)
                
            one_d_i = [x.reshape(-1) for x in array_of_indices]
            for l, x in zip(one_dimensional_array_of_indices, one_d_i):
                l.append(x)
        
        one_dimensional_values = [concatenate(x) for x in one_dimensional_values]
        one_dimensional_array_of_indices = [concatenate(x) for x in one_dimensional_array_of_indices]
        
        for x in extra_attributes:
            extra_kw[x] = concatenate(extra_kw[x])
        selected_setters = list([setter for setter in atss[0].select_setters_for(original_attributes)])
        def next_request(index, setters):
            if index < len(setters):
                setter = setters[index]
                return setter.set_attribute_values_async(atss[0], attributes, one_dimensional_values, *one_dimensional_array_of_indices, **extra_kw)
            else:
                return None
        
        request = ASyncRequestSequence(next_request, args = (selected_setters,))
        return request

class EvolveHydrodynamicsCodeWithPeriodicBoundaries(object):
    def __init__(self, code):
        self.code = code
    
    def set_parameters(self):
        self.code.parameters.x_boundary_conditions = ("periodic","periodic")
        self.code.parameters.y_boundary_conditions = ("periodic","periodic")
        self.code.stopping_conditions.number_of_steps_detection.disable()
        
    def init_channels(self):
        pass 
    
    def evolve_model(self, time, endtime):
        while self.code.model_time < time:
                
            if (self.code.get_timestep() + self.code.model_time) >= time and time == endtime:
                self.code.parameters.must_evolve_to_exact_time = True
            
            self.code.evolve_model(time)


class CalculateLinearWave1D(object):
    gamma = 5.0/3.0
    wave_flag = 0
    
    def __init__(self, 
        number_of_grid_points =  10, 
        number_of_workers = 1, 
        name_of_the_code = "athena",
        amplitude = 1e-6 | speed,
        vflow_factor = 1.0,
        grid_length = 1.0 | length,
        number_of_steps = 10,
        use_boundaries = True,
        number_of_codes = 2
    ):
        
        self.number_of_grid_points = number_of_grid_points
        self.number_of_workers = number_of_workers
        self.name_of_the_code = name_of_the_code
        self.amplitude = amplitude
        self.vflow_factor = vflow_factor
        self.grid_length = 1.0 | length
        self.number_of_steps = number_of_steps
        self.number_of_codes = number_of_codes
        self.dimensions_of_mesh = (
            self.number_of_grid_points, 
            self.number_of_grid_points, 
            1
        )
        self.nghost = 4
        self.use_boundaries = use_boundaries
        
    def new_instance_of_code(self):
        if 0:
            if self.name_of_the_code == 'athena':
                self.name_of_the_code = 'capreole'
            else:
                self.name_of_the_code = 'athena'
            
        attribute = "new_instance_of_{0}_code".format(self.name_of_the_code.lower())
        return getattr(self,attribute)()
        
    def new_instance_of_capreole_code(self):
        result=Capreole(number_of_workers=self.number_of_workers)
        self.dimensions_of_mesh =  (
            self.number_of_grid_points, 
            self.number_of_grid_points, 
            4
        )
        self.nghost = 2
        return result
        
    def new_instance_of_athena_code(self):
        from amuse.community.athena.interface import Athena
        result=Athena(number_of_workers=self.number_of_workers)
        result.parameters.gamma = self.gamma
        result.parameters.courant_number=0.4
        self.nghost = 4
        return result
        

    def new_instance_of_mpiamrvac_code(self):
        raise Exception("MPIAMRVAC does not yet have support for detailed boundaries in amuse")
        from amuse.community.mpiamrvac.interface import MpiAmrVac
        result=MpiAmrVac(mode="2d", number_of_workers=self.number_of_workers, debugger="xterm")
        result.set_parameters_filename(result.default_parameters_filename)
        result.initialize_code()
        return result
        
    def set_parameters(self, codes, is_live_code, evolve):
        
        mesh_for_code = list(self.dimensions_of_mesh)
        i = 1
        for is_live, instance in zip(is_live_code, codes):
            if not is_live:
                continue
            local_mesh = list(mesh_for_code)
            #local_mesh[1] *= i #(1 + (i%2))
            instance.parameters.mesh_size = list(local_mesh)
            i += 1
            
            instance.parameters.length_x = self.grid_length
            instance.parameters.length_y =  self.grid_length / self.number_of_codes
            instance.parameters.length_z = self.grid_length
            
            instance.parameters.z_boundary_conditions = ("periodic","periodic")
            
        evolve.set_parameters()
        
        for is_live, instance in zip(is_live_code, codes):
            if not is_live:
                continue
            instance.commit_parameters()
    
    def new_grid(self):
        grid = Grid.create(self.dimensions_of_mesh, [1,1,1] | length)
        self.clear_grid(grid)
        return grid
    
    def clear_grid(self, grid):
        density = mass / length**3
        momentum =  speed * density
        energy =  mass / (time**2 * length)

        grid.rho =  0.0 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.0 | energy
    
        return grid
    
    
    def new_rhoe_right_eigenmatrix(self, velocity, amplitude, enthalpy):
        right_eigenmatrix = numpy.zeros ( (5,5) ) | speed
        
        right_eigenmatrix[0][0] = 1.0 | speed;
        right_eigenmatrix[1][0] = velocity[0] - amplitude;
        right_eigenmatrix[2][0] = velocity[1];
        right_eigenmatrix[3][0] = velocity[2];
        right_eigenmatrix[4][0] = (1.0 | time/length) * (enthalpy - velocity[0]*amplitude);
        #right_eigenmatrix[0][1] = 0.0;
        #right_eigenmatrix[1][1] = 0.0;
        right_eigenmatrix[2][1] = 1.0 | speed;
        #right_eigenmatrix[3][1] = 0.0; 
        right_eigenmatrix[4][1] = velocity[1];

        #right_eigenmatrix[0][2] = 0.0; */
        #right_eigenmatrix[1][2] = 0.0; */
        #right_eigenmatrix[2][2] = 0.0; */
        right_eigenmatrix[3][2] = 1.0 | speed;
        right_eigenmatrix[4][2] = velocity[2];

        right_eigenmatrix[0][3] = 1.0 | speed;
        right_eigenmatrix[1][3] = velocity[0];
        right_eigenmatrix[2][3] = velocity[1];
        right_eigenmatrix[3][3] = velocity[2];
        right_eigenmatrix[4][3] = 0.5*velocity.length();

        right_eigenmatrix[0][4] = 1.0 | speed;
        right_eigenmatrix[1][4] = velocity[0] + amplitude;
        right_eigenmatrix[2][4] = velocity[1];
        right_eigenmatrix[3][4] = velocity[2];
        right_eigenmatrix[4][4] = (1.0 | time/length) *  (enthalpy + velocity[0]*amplitude);
        return right_eigenmatrix

    def initialize_grid(self, grid):
        density = mass / length**3
        momentum =  speed * density
        energy =  mass / (time**2 * length)
        rho =  1.0 | density
        pressure = (1.0/self.gamma) | (mass / (length * time**2))
        vx = (self.gamma * pressure / rho).sqrt()
        velocity = self.vflow_factor * vx * [1.0, 0.0, 0.0]
        velocity_squared = velocity.length()
        energy = (pressure/(self.gamma - 1.0) + (0.5 | length / time )*rho*velocity_squared)
        enthalpy = (energy + pressure)/rho;
        amplitude_squared = (self.gamma - 1.0) * max(enthalpy - (0.5 | length/time)* velocity_squared, 1e-100 | enthalpy.unit) 
        amplitude =amplitude_squared.sqrt()
        
        nwave = 5
        eigenvalues = ([0] * nwave) | speed
        eigenvalues[0] = velocity[0] - amplitude
        eigenvalues[1] = velocity[0]
        eigenvalues[2] = velocity[0]
        eigenvalues[3] = velocity[0]
        eigenvalues[4] = velocity[0] + amplitude
        
        right_eigenmatrix = self.new_rhoe_right_eigenmatrix(velocity, amplitude, enthalpy)
        
        grid.rho = rho
        grid.energy = energy
        grid.rhovy = rho*self.vflow_factor*(1.0 | speed)
        
        wave = self.amplitude*numpy.sin(grid.y * (2.0 | length**-1)*numpy.pi)
        grid.rho += wave*right_eigenmatrix[0][self.wave_flag] * (1.0 |mass * time**2 / length**5)
        grid.rhovx += wave*right_eigenmatrix[3][self.wave_flag] * (1.0 |mass * time / length**4)
        grid.rhovy += wave*right_eigenmatrix[1][self.wave_flag] * (1.0 |mass * time / length**4)
        grid.rhovz += wave*right_eigenmatrix[2][self.wave_flag] * (1.0 |mass * time / length**4)
        grid.energy += wave*right_eigenmatrix[4][self.wave_flag] *(1.0 | mass  / length**3)
        
    def store_grids(self, grids, step):
        if __name__ == '__plot__':
            return
        
        grids_in_memory = [x.copy() for x in grids]
        io.write_set_to_file(
            grids_in_memory, 
            "linear_wave_{2}_{0}_{1}.vtu".format(self.number_of_grid_points, step, self.name_of_the_code),
            "vtu",
            is_multiple=True
        )
    
    def new_evolve_object(self, instances, offsets, is_live):
        """Returns a special object to evolve to code in time"""
        if self.use_boundaries:
            #return EvolveHydrodynamicsCodeWithAmusePeriodicBoundaries(instance[0], self.number_of_grid_points, self.nghost)
            #return EvolveHydrodynamicsCodeWithAmusePeriodicBoundariesAndNCodes(instance, self.number_of_grid_points, self.number_of_grid_points / self.number_of_codes , self.nghost)
            nodes = []
            for offset, code, live in zip(offsets,instances,is_live):
                node = Node(code, offset, not live)
                nodes.append(node)
            return EvolveHydrodynamicsCodeWithAmusePeriodicBoundariesAndNCodesWithDifferentGridsSizes(nodes)
        else:
            return EvolveHydrodynamicsCodeWithPeriodicBoundaries(instance[0])
            
    def initialize(self):
        self.all_codes = []
        self.codes = []
        self.offsets = []
        self.is_live = []
        
        y = 0 *  self.grid_length
        dy = [0,1,0] * (self.grid_length / self.number_of_codes)
        dx = [1,0,0] * (self.grid_length)
        offset = [0,0,0] * self.grid_length
        for i in range(self.number_of_codes):
            code = self.new_instance_of_code()
            print i, id(code)
            self.codes.append(code)
            self.offsets.append(offset)
            self.is_live.append(True)
            offset += dy
        
        self.offsets.append(self.offsets[0] - dy)
        self.codes.append(self.codes[self.number_of_codes-1])
        self.is_live.append(False)
        
        self.offsets.append(self.offsets[self.number_of_codes-1] + dy)
        self.codes.append(self.codes[0])
        self.is_live.append(False)
        
        #periodical in x, one code length
        for i in range(self.number_of_codes+2):
            code = self.codes[i]
            offset = self.offsets[i]
            self.offsets.append(offset - dx)
            self.codes.append(code)
            self.is_live.append(False)
            
        #periodical in x, one code length
        for i in range(self.number_of_codes+2):
            code = self.codes[i]
            offset = self.offsets[i]
            self.offsets.append(offset + dx)
            self.codes.append(code)
            self.is_live.append(False)
        
        
            
        self.evolve = self.new_evolve_object(self.codes, self.offsets, self.is_live)
        
        self.set_parameters(self.codes, self.is_live, self.evolve)
        
        self.start_grids = []
        
        offset = 0.0 * self.grid_length
        for is_live, instance in zip(self.is_live, self.codes):
            if not is_live:
                continue
            for x in instance.itergrids():
                inmem = x.copy()
                self.clear_grid(inmem)
                inmem.y += offset
                self.initialize_grid(inmem)
                self.start_grids.append(inmem)
                from_model_to_code = inmem.new_channel_to(x)
                from_model_to_code.copy()
            offset += self.grid_length / self.number_of_codes
            
        self.get_grids()
        for is_live, x in zip(self.is_live, self.codes):
            if not is_live:
                continue
            x.initialize_grid()
        self.evolve.init_channels()
        
    def evolve_model(self,time, endtime):
        for code in self.codes:
            code.parameters.must_evolve_to_exact_time = False
        self.evolve.evolve_model(time, endtime)
            
    
    def get_grids(self):
        result = []
        
        offset = 0.0 * self.grid_length
        for is_live, code in zip(self.is_live, self.codes):
            if not is_live:
                continue
            for x in code.itergrids():
                inmem = x.copy()
                inmem.y += offset
                result.append(inmem)
            offset += self.grid_length / self.number_of_codes

        return result
    
    def stop(self):
        print "terminating code"
        
        for code in self.codes:
            code.stop()

import sys
def main():
    number_of_grid_points = 40
    name_of_the_code = 'athena'
    number_of_steps = 2000
    vflow_factor = -1.0
    pertubation_amplitude = 1e-4 | speed
    grid_length = 1.0 | length
    number_of_codes = int(sys.argv[1])
    #if number_of_grid_points % number_of_codes != 0:
    #    raise Exception("grid points should be dividable by the number of codes")
    model1 = CalculateLinearWave1D(
        number_of_grid_points = number_of_grid_points,
        number_of_workers = 1,
        name_of_the_code = name_of_the_code,
        amplitude = pertubation_amplitude,
        vflow_factor = vflow_factor,
        grid_length = grid_length,
        number_of_steps = number_of_steps,
        use_boundaries = True,
        number_of_codes = number_of_codes
    )
    if 0:
        model2 = CalculateLinearWave1D(
            number_of_grid_points = number_of_grid_points,
            number_of_workers = 1,
            name_of_the_code = name_of_the_code,
            amplitude = pertubation_amplitude,
            vflow_factor = vflow_factor,
            grid_length = grid_length,
            number_of_steps = number_of_steps,
            use_boundaries = False,
            number_of_codes = number_of_codes
        )
        
    if not IS_PLOT_AVAILABLE:
        print "Plot is not available. stop now"
        return
    model1.initialize()
    #model2.initialize()
    
    
    grids1 = model1.get_grids()
    #grids2 = model2.get_grids()

    
    figure = pyplot.figure(figsize=(10,5))
    plot1 = figure.add_subplot(1,1,1)
    lines = []
    ys = []
    colors = ['b','g','r','c','m']
    for i,grid in enumerate(grids1):
        y = grid.y[0,...,0].value_in(length)
        ys.append(y)
        rho = grid.rho[0,...,0].value_in(density)
        print rho-1.0,colors[i % len(colors)]
        line = plot1.scatter(y,rho, c = colors[i % len(colors)])#[0]
        lines.append(line)
    code0 = model1.codes[0]
    boundary0 = code0.get_boundary_grid("ybound2")
    y = boundary0.y[0,...,0].value_in(length)
    rho = boundary0.rho[0,...,0].value_in(density)
    lines.append(plot1.scatter(y,rho, c = 'y'))
    
    plot1.set_xlim(0,1)
    plot1.set_ylim(1-pertubation_amplitude.value_in(speed), 1+pertubation_amplitude.value_in(speed))
    end_time = 10.0 | time
    dt = end_time / number_of_steps
    
    t = dt
    step = 1
    
    title = figure.suptitle('{0:.3f}'.format(0.0))
    variables = [t, step]

    def update(xx):
        t, step = variables
        #if step >= 1: return
        title.set_text('{0:.3f}'.format(t.value_in(time)))
        model1.evolve_model(t, end_time)
        #model2.evolve_model(t, end_time)
        t += dt
        grids1 = model1.get_grids()
        #grids2 = model2.get_grids()
        
        for line, grid in zip(lines, grids1):
            y = grid.y[0,...,0].value_in(length)
            rho = grid.rho[0,...,0].value_in(density)
            offsets = numpy.vstack([y,rho]).transpose()
            
            #print line.get_offsets()
            line.set_offsets(offsets)
            #line.set_data(y,rho)
            #line = plot1.plot(y,rho)[0]
            #lines.append(line)
        
        code0 = model1.codes[0]
        boundary0 = code0.get_boundary_grid("ybound2")
        y = boundary0.y[0,...,0].value_in(length)
        rho = boundary0.rho[0,...,0].value_in(density)
        offsets = numpy.vstack([y,rho]).transpose()
        lines[-1].set_offsets(offsets)
        print t
        step += 1
        variables[0] = t
        variables[1] = step
        return lines
    if 1:
        process = animation.FuncAnimation(
            figure, 
            update, 
            numpy.arange(1, 200), 
            interval=2000,
            blit=False
        )
    else:
        update(0)
        #update(0)
        pass
    
        
    pyplot.show()
    model1.stop()
    #model2.stop()
    
if __name__ == "__main__":
    main()
