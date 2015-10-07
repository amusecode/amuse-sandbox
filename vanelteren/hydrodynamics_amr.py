"""
In this script we simulate a 1d linear wave in a 2d field, 
the periodic boundaries are not handled in the code but by amuse
(much slower at time of writing)
"""
import numpy
import math
import time as T
from itertools import groupby

from matplotlib import pyplot
from mpl_toolkits.mplot3d import axes3d

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

from amuse.datamodel import SubGrid
from amuse.datamodel import Grid
from amuse.datamodel import GridPoint
from amuse.datamodel import AbstractGrid
from amuse.datamodel import indexing
from amuse.datamodel import grid_attributes

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
    
def amr_wave_filter(level):
    return 1.0e-2
    
def lohner_error_estimation(grid, attributes):
    result = 0.0
    shape = grid.shape
    dimensions = 1
    if shape[1] > 1:
         dimensions += 1
    if shape[2] > 1:
         dimensions += 1
         
    for attribute in attributes:
        value = getattr(grid, attribute)
        numerator=None
        denominator=None
        for dimension in range(dimensions):
            imin0 = 0 if dimension == 0 else 1
            imax0 = -2 if dimension == 0 else -1
            imin1 = 2 if dimension == 0 else 1
            imax1 = None if dimension == 0 else -1
            
            if dimensions > 1:
                jmin0 = 0 if dimension == 1 else 1
                jmax0 = -2 if dimension == 1 else -1
                jmin1 = 2 if dimension == 1 else 1
                jmax1 = None if dimension == 1 else -1
            else:
                jmin0 = jmax0 = jmin1 = jmax1 = None
                
            if dimensions > 2:
                kmin0 = 0 if dimension == 2 else 1
                kmax0 = -2 if dimension == 2 else -1
                kmin1 = 2 if dimension == 2 else 1
                kmax1 = None if dimension == 2 else -1
            else:
                kmin0 = kmax0 = kmin1 = kmax1 = None
            
            dvalue = value[imin1:imax1,jmin1:jmax1,kmin1:kmax1] - value[imin0:imax0,jmin0:jmax0,kmin0:kmax0]
            
            for dimension2 in range(dimensions):
                
                imin02 = 0 if dimension2 == 0 else 1
                imax02 = -2 if dimension2 == 0 else -1
                imin12 = 2 if dimension2 == 0 else 1
                imax12 = None if dimension2 == 0 else -1
                
                if dimensions > 1:
                    jmin02 = 0 if dimension2 == 1 else 1
                    jmax02 = -2 if dimension2 == 1 else -1
                    jmin12 = 2 if dimension2 == 1 else 1
                    jmax12 = None if dimension2 == 1 else -1
                else:
                    jmin02 = jmax02 = jmin12 = jmax12 = None
                
                if dimensions > 2:
                    kmax02 = -2 if dimension2 == 2 else -1
                    kmin02 = 0 if dimension2 == 2 else 1
                    kmin12 = 2 if dimension2 == 2 else 1
                    kmax12 = None if dimension2 == 2 else -1
                else:
                    kmin02 = kmax02 = kmin12 = kmax12 = None
                
                ddvalue = dvalue[imin12:imax12,jmin12:jmax12,kmin12:kmax12] - dvalue[imin02:imax02,jmin02:jmax02,kmin02:kmax02]
                if numerator is None:
                    numerator = (ddvalue * ddvalue)
                else:
                    numerator += (ddvalue * ddvalue)
                
        for dimension in range(dimensions):
            absvalue = abs(value)
            
            imin0 = 0 if dimension == 0 else 1
            imax0 = -2 if dimension == 0 else -1
            imin1 = 2 if dimension == 0 else 1
            imax1 = None if dimension == 0 else -1
            
            if dimensions > 1:
                jmin0 = 0 if dimension == 1 else 1
                jmax0 = -2 if dimension == 1 else -1
                jmin1 = 2 if dimension == 1 else 1
                jmax1 = None if dimension == 1 else -1
            else:
                jmin0 = jmax0 = jmin1 = jmax1 = None
            
            if dimensions > 2:
                kmin0 = 0 if dimension == 2 else 1
                kmax0 = -2 if dimension == 2 else -1
                kmin1 = 2 if dimension == 2 else 1
                kmax1 = None if dimension == 2 else -1
            else:
                kmin0 = kmax0 = kmin1 = kmax1 = None
                
            sum_value= absvalue[imin1:imax1,jmin1:jmax1,kmin1:kmax1] + absvalue[imin0:imax0,jmin0:jmax0,kmin0:kmax0]
            
            imin0 = 0 if dimension == 0 else 2
            imax0 = -4 if dimension == 0 else -2
            imin1 = 4 if dimension == 0 else 2
            imax1 = None if dimension == 0 else -2
            
            if dimensions > 1:
                jmin0 = 0 if dimension == 1 else 2
                jmax0 = -4 if dimension == 1 else -2
                jmin1 = 4 if dimension == 1 else 2
                jmax1 = None if dimension == 1 else -2
            else:
                jmin0 = jmax0 = jmin1 = jmax1 = None
            
            if dimensions > 2:
                kmin0 = 0 if dimension == 2 else 2
                kmax0 = -4 if dimension == 2 else -2
                kmin1 = 4 if dimension == 2 else 2
                kmax1 = None if dimension == 2 else -2
            else:
                kmin0 = kmax0 = kmin1 = kmax1 = None
            
            if dimensions == 1:
                field_value = value[2:-2,:,:]
            elif dimensions == 2:
                field_value = value[2:-2,2:-2,:]
            elif dimensions == 3:
                field_value = value[2:-2,2:-2,2:-2]
            
            diff_value = abs(value[imin1:imax1,jmin1:jmax1,kmin1:kmax1] - field_value) +  abs(field_value - value[imin0:imax0,jmin0:jmax0,kmin0:kmax0])
            
            
            for dimension2 in range(dimensions):
                
                imin02 = 0 if dimension2 == 0 else 1
                imax02 = -2 if dimension2 == 0 else -1
                imin12 = 2 if dimension2 == 0 else 1
                imax12 = None if dimension2 == 0 else -1
                
                    
                if dimensions > 1:
                    jmin02 = 0 if dimension2 == 1 else 1
                    jmax02 = -2 if dimension2 == 1 else -1
                    jmin12 = 2 if dimension2 == 1 else 1
                    jmax12 = None if dimension2 == 1 else -1
                else:
                    jmin02 = jmax02 = jmin12 = jmax12 = None
                
                if dimensions > 2:
                    kmin02 = 0 if dimension2 == 2 else 1
                    kmax02 = -2 if dimension2 == 2 else -1
                    kmin12 = 2 if dimension2 == 2 else 1
                    kmax12 = None if dimension2 == 2 else -1
                else:
                    kmin02 = kmax02 = kmin12 = kmax12 = None
                 
                ssvalue = sum_value[imin12:imax12,jmin12:jmax12,kmin12:kmax12]  + sum_value[imin02:imax02,jmin02:jmax02,kmin02:kmax02]
                tmp = (diff_value + (amr_wave_filter(1) * ssvalue))**2
                if denominator is None:
                    denominator = tmp
                else:
                    denominator += tmp
        result += numpy.sqrt(numerator / (denominator.maximum(1e-6 | denominator.unit)))
    return result
        
class LohnenErrorEstimationTests(TestCase):
    def test1(self):
        input = Grid.create((7,6,5), [7,6,5] | generic_unit_system.length)
        input.rho = 0.0 | generic_unit_system.density
        error = lohner_error_estimation(input, ["rho",])
        self.assertEquals(error, 0)
        
    def test2(self):
        from amuse.community.mpiamrvac.interface import MpiAmrVac
        code = MpiAmrVac(mode="1d", redirection = "none")
        code.parameters.x_boundary_conditions = ("periodic","periodic")
        code.parameters.mesh_length = (12.0,1,1) | generic_unit_system.length
        code.parameters.mesh_size = (40, 1, 1)
        code.parameters.maximum_number_of_grid_levels = 2
        for grid in code.itergrids():
            copy = grid.copy()
            copy.rho = (1.0 + (0.5 * (numpy.sin(copy.x/(12.0 | generic_unit_system.length)*  2 * numpy.pi)))) | generic_unit_system.density
            channel = copy.new_channel_to(grid)
            channel.copy()
        must_refine = code.refine_grid()
        print must_refine
            
            
        grid = Grid.create((84,1,1), [12,12,12] | generic_unit_system.length)
        r = ((grid.x**2  + grid.y**2)).sqrt()
        grid.rho = 1.0 + 0.5 * (numpy.sin(grid.x/(12.0 | generic_unit_system.length)*  2 * numpy.pi)) | generic_unit_system.density
        error = lohner_error_estimation(grid, ["rho",])
        
        rho = grid[2:-2,:,:].rho[...,...,0].value_in(generic_unit_system.density)
        rho = error[...,...,0]
        print rho
        figure = pyplot.figure(figsize=(6,6))
        plot = figure.add_subplot(1,1,1)#, projection='3d')
        plot.plot(
            grid[2:-2,:,:].x[...,...,0].value_in(generic_unit_system.length),
            rho
        )
        figure.savefig('test2.png')
        pyplot.show()
        
        #self.assertEquals(error, 0)
        

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

class OffsetGrid(AbstractGrid):
    DEFAULT_AXES_NAMES = ('x', 'y', 'z')
    
    
    def __init__(self, grid, offset, **kwargs):
        AbstractGrid.__init__(self)
        self._private.offset = offset
        self._private.grid = grid
        if "axes_names" in kwargs and (not kwargs['axes_names'] is None):
            self._private.axes_names = kwargs['axes_names']
        else:
            self._private.axes_names = self.DEFAULT_AXES_NAMES
        self._private.axes_indices = {}
        for i, axis in enumerate(self._private.axes_names):
            self._private.axes_indices[axis] = i
            
        self.add_vector_attribute("position", self._private.axes_names[0:len(self._private.grid.shape)])
        
    def can_extend_attributes(self):
        return self._private.grid.can_extend_attributes()
            
    def get_values_in_store(self, indices, attributes, by_key = True):
        result = self._private.grid.get_values_in_store(indices, attributes)
        
        offset_result = []
        for x, value in zip(attributes, result):
            if x in self._private.axes_indices:
                value += self._private.offset[self._private.axes_indices[x]]
            offset_result.append(value)
            
        return offset_result
        
    def set_values_in_store(self, indices, attributes, values, by_key = True):
        
        offset_values = []
        for x, value in zip(attributes, values):
            if x in self._private.axes_indices:
                value -= self._private.offset[self._private.axes_indices[x]]
            offset_values.append(value)
            
        self._private.grid.set_values_in_store(indices, attributes, offset_values)
        
    def set_values_in_store_async(self, indices, attributes, values, by_key = True):
        
        offset_values = []
        for x, value in zip(attributes, values):
            if x in self._private.axes_indices:
                value -= self._private.offset[self._private.axes_indices[x]]
            offset_values.append(value)
        
        return self._private.grid.set_values_in_store_async(indices, attributes, offset_values)
        

    def get_attribute_names_defined_in_store(self):
        return self._private.grid.get_attribute_names_defined_in_store()
        
    def get_defined_settable_attribute_names(self):
        return self._private.grid.get_defined_settable_attribute_names()

    def _original_set(self):
        return self
        
    def get_all_keys_in_store(self):
        return self._private.grid.get_all_keys_in_store()
    
    def __getitem__(self, index):
        if indexing.number_of_dimensions_after_index(self.number_of_dimensions(), index) == 0:
            return GridPoint(index, self)
        else:
            return SubGrid(self, index)
    
    def iter_cells(self):
        shape = numpy.asarray(self.shape)
        
        index = 0 * shape
        
        while index[0] < shape[0]:
            yield self._get_gridpoint(tuple(index))
            
            index[-1] += 1
            for i in range(len(self.shape) - 1, 0, -1):
                if index[i] >= shape[i]:
                    index[i] = 0
                    index[i-1] += 1
            
                
    def _get_gridpoint(self, index):
        return GridPoint(index, self)
        
    def number_of_dimensions(self):
        return len(self.shape)
        
    @property
    def shape(self):
        return self._private.grid.shape
        
    @property
    def size(self):
        return numpy.prod(self.shape)
        
    def indices(self):
        return numpy.indices(self.shape)
        

class OffsetGridTests(TestCase):
    def test1(self):
        input = Grid.create((5,4,3), [5,4,3] | generic_unit_system.length)
        self.assertEquals(input[0][0][0].x , 0.5 | generic_unit_system.length)
        self.assertEquals(input[0][0][0].position , (0.5,0.5,0.5) | generic_unit_system.length)
        offset = OffsetGrid(input, (2, 3, 4) | generic_unit_system.length)
        self.assertEquals(offset[0][0][0].x,  2.5 | generic_unit_system.length)
        self.assertEquals(offset[0][0][0].position , (2.5,3.5,4.5) | generic_unit_system.length)
        self.assertEquals(offset[0:5,0,0].position , [(2.5,3.5,4.5),(3.5,3.5,4.5),(4.5,3.5,4.5),(5.5,3.5,4.5),(6.5,3.5,4.5)] | generic_unit_system.length)
        
        
    

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
        
class HydrodynamcisAMR(object):
    
    def __init__(self, mesh_size, mesh_length, create_new_code = None):
        self.create_new_code = create_new_code
        self.mesh_size = mesh_size
        self.mesh_length = mesh_length * 1.0
        self.cellsize = self.mesh_length / self.mesh_size
        self.min_timestep = None
        self.channels = []
        self.nodes = []
        self.max_levels = 4
        
        self.step = 0
    
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
                    channel = HydroStateChannel(other_code, node.code, boundary[selection], offset - other_offset, factor)
                    channels.append(channel)
            
                        
        self.channels = list(reversed(sorted(channels, key = lambda x : x.factor)))
        
        self.channels_grouped_by_factor_and_code = []
        self.channels_grouped_by_factor_and_to_code = []
        for factor, group in groupby(self.channels, lambda x : x.factor):
            grouped_by_factor = list(group)
            grouped_by_factor = list(sorted(grouped_by_factor, key = lambda x : id(x.from_code)))
            for from_code_id, code_group in groupby(grouped_by_factor, lambda x : id(x.from_code)): 
                grouped_by_factor_and_code = list(code_group)
                self.channels_grouped_by_factor_and_code.append(grouped_by_factor_and_code)
                
                
            grouped_by_factor = list(sorted(grouped_by_factor, key = lambda x : id(x.to_code)))
            for to_code_id, code_group in groupby(grouped_by_factor, lambda x : id(x.to_code)): 
                grouped_by_factor_and_to_code = list(code_group)
                grouped_by_factor_and_to_code = list(sorted(grouped_by_factor_and_to_code, key = lambda x : x.copy_grid.size))
                self.channels_grouped_by_factor_and_to_code.append(grouped_by_factor_and_to_code)
            
        self.copy_to_boundary_cells()
        
    def copy_to_boundary_cells(self):
        print "copy to boundary..."
        if 0:
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
            #print len(self.channels)
            for channel in self.channels:
                channel.copy()
        
        print "copy to boundary...done"
    
    def setup_codes(self):
        
        # start with one code covering the entire mesh
        mesh_size = numpy.asarray(self.mesh_size)
        mesh_length = self.mesh_length
        block_size = 20
        overflow_per_dimension = mesh_size % block_size
        print overflow_per_dimension, mesh_size, block_size
        for x,y in zip(mesh_size, overflow_per_dimension):
            if x>1 and y != 0:
                raise Exception("cannot divide the cells in dimension")
        number_of_blocks_per_dimension = numpy.maximum(mesh_size / block_size, (1,1,1))
        
        deltax = mesh_length / number_of_blocks_per_dimension
        node_size = numpy.asarray((block_size, 1, 1))
        if self.mesh_size[1] > 1:
            node_size[1] = block_size
        if self.mesh_size[2] > 1:
            node_size[2] = block_size
        print node_size
        #define the nodes
        self.nodes = []
        
        for i in range(number_of_blocks_per_dimension[0]):
            offset_x = i * deltax[0]
            for j in range(number_of_blocks_per_dimension[1]):
                offset_y = j * deltax[1]
                for k in range(number_of_blocks_per_dimension[2]):
                    offset_z = k * deltax[2]
                    code = self.create_new_code()
                    code.parameters.mesh_size = node_size
                    code.parameters.mesh_length = deltax
                    code.parameters.x_boundary_conditions = ("interface","interface")
                    if self.mesh_size[1] > 1:
                        code.parameters.y_boundary_conditions = ("interface","interface")
                    else:
                        code.parameters.y_boundary_conditions = ("periodic","periodic")
                    if self.mesh_size[2] > 1:
                        code.parameters.z_boundary_conditions = ("interface","interface")
                    else:
                        code.parameters.z_boundary_conditions = ("periodic","periodic")
                    code.stopping_conditions.number_of_steps_detection.enable()
        
                    code.commit_parameters()
                    offset = deltax * (i,j,k)
                    print deltax, offset
                    node = Node(code, offset, False)
                    node.setup()
                    self.nodes.append(node)
        print len(self.nodes)
        ghost_nodes = self.new_ghost_nodes(self.nodes)
        self.nodes.extend(ghost_nodes)
            
            
    def evolve_model(self, time, endtime = None):
        if endtime is None:
            endtime = time
            
        code = self.nodes[0].code
        t_unit = code.get_timestep().unit
        while ((code.model_time - time)/time) < -1e-14:
            timestep = self.get_timestep_from_codes()
            self.set_timestep(code.model_time, time, endtime, timestep)
           
            self.evolve_all_codes(time)
            
            self.copy_to_boundary_cells()
            
            self.step += 1
            if self.step % 5 == 0:
                print "refining the grid...", self.step,  code.model_time
                has_refined = self.refine_grid()
                print "...done refined = ", has_refined
                code = self.nodes[0].code
    
    def get_timestep_from_codes(self):
        timesteps = [] | t_unit
            for x in self.livenodes():
                timesteps.append(x.code.get_timestep())
        return timesteps.min()
        
    def set_timestep(self, current_time, nexttime, endtime, timestep):
        if (timestep + current_time) >= nexttime and nexttime == endtime:
            for x in self.livenodes():
                x.code.parameters.must_evolve_to_exact_time = True
        for x in self.livenodes():
            x.code.set_timestep(timestep)
    
    def evolve_all_codes(self):
        pool = AsyncRequestsPool()
        for x in self.livenodes():
            request = x.code.evolve_model.async(time)
            pool.add_request(request, lambda x : x.result())
        pool.waitall()
            
    def livenodes(self):
        for x in self.nodes:
            if not x.is_ghost:
                yield x
                
    def itergrids(self):
        for x in self.livenodes():
            for grid in x.code.itergrids():
                yield OffsetGrid(grid, x.offset)
                    
    def refine_grid(self):
        refined_nodes = []
        is_refined = False
        for x in self.nodes:
            if x.is_ghost:
                continue
            can_refine = numpy.all(x.cellsize > (self.cellsize / 2**self.max_levels))
            if not can_refine: 
                refined_nodes.append(x)
                continue
            must_refine = False
            for grid in x.code.itergrids():
                copy = grid.copy()
                error = lohner_error_estimation(copy, ["rho", "rhovx", "rhovy", "rhovz", "energy"])
                if numpy.any(error > 0.5):
                    must_refine |= True
            if must_refine:
                if 1:
                    is_refined = True
                    factor = (2.,1.,1.)
                    maxi = maxj = maxk = 1
                    maxi = 2
                    if self.mesh_size[1] > 1:
                        factor = (2.,2.,1.)
                        maxj = 2
                    if self.mesh_size[2] > 1:
                        factor = (2.,2.,2.)
                        maxk = 2
                    print factor
                    
                    for i in range(maxi):
                        for j in range(maxj):
                            for k in range(maxk):
                                code = self.create_new_code()
                                code.parameters.mesh_size = numpy.asarray(self.mesh_size)
                                code.parameters.mesh_length = x.size / factor
                                code.parameters.x_boundary_conditions = ("interface","interface")
                                if x.number_of_grid_points[1] > 1:
                                    code.parameters.y_boundary_conditions = ("interface","interface")
                                else:
                                    code.parameters.y_boundary_conditions = ("periodic","periodic")
                                if x.number_of_grid_points[2] > 1:
                                    code.parameters.z_boundary_conditions = ("interface","interface")
                                else:
                                    code.parameters.z_boundary_conditions = ("periodic","periodic")
                                code.stopping_conditions.number_of_steps_detection.enable()
                    
                                code.commit_parameters()
                                for grid in code.itergrids():
                                    channel = HydroStateChannel(x.code, code, grid, (x.size * (i*0.5, j*0.5, k*0.5)), 1)
                                    channel.copy()
                                
                                node = Node(code, x.offset + (x.size * (i*0.5, j*0.5, k*0.5)), False)
                                node.setup()
                                refined_nodes.append(node)
                    x.code.stop()
                else:
                    refined_nodes.append(x)
            else:
                refined_nodes.append(x)
    
        if is_refined:
            ghost_nodes = self.new_ghost_nodes(refined_nodes)
            nodes = []
            nodes.extend(refined_nodes)
            nodes.extend(ghost_nodes)
            print len(nodes)
            self.nodes = nodes
            self.init_channels()
        return is_refined
    
    def new_ghost_nodes(self, refined_nodes):
        mesh_length = self.mesh_length
        ghost_nodes = []
        for x in refined_nodes:
            code = x.code
            node_size = x.size
            node_offset = x.offset
            print node_offset
            for i in range(2):
                print i, node_offset[i] == quantities.zero, node_offset[i] + node_size[i] == mesh_length[i]
                
            
            if node_offset[0] == quantities.zero:
                ghost_nodes.append(Node(code, mesh_length * (1, 0, 0) + node_offset * (0,1,1), True )) 
            if node_offset[0] + node_size[0] == mesh_length[0]:
                ghost_nodes.append(Node(code, node_size * (-1, 0, 0) + node_offset * (0,1,1) , True )) 
            
            if x.number_of_grid_points[1] > 1:
                if node_offset[1] == quantities.zero:
                    ghost_nodes.append(Node(code, mesh_length * (0, 1, 0) + node_offset * (1,0,1), True )) 
                if node_offset[1] + node_size[1] == mesh_length[1]:
                    ghost_nodes.append(Node(code, node_size * (0, -1, 0)  + node_offset * (1,0,1), True ))  
                if node_offset[0] == quantities.zero and node_offset[1] == quantities.zero:
                    ghost_nodes.append(Node(code, mesh_length * (1, 1, 0) + node_offset * (0,0,1), True )) 
                if node_offset[0] == quantities.zero and node_offset[1] + node_size[1] == mesh_length[1]:
                    ghost_nodes.append(Node(code, mesh_length * (1, 0, 0) + (node_size * (0, -1, 0)) + node_offset * (0,0,1), True )) 
                if node_offset[1] == quantities.zero and node_offset[0] + node_size[0] == mesh_length[0]:
                    ghost_nodes.append(Node(code, mesh_length * (0, 1, 0) + (node_size * (-1, 0, 0)) + node_offset * (0,0,1), True )) 
                if node_offset[1] + node_size[1] == mesh_length[1] and node_offset[0] + node_size[0] == mesh_length[0]:
                    ghost_nodes.append(Node(code, node_size * (-1, -1, 0) + node_offset * (0,0,1), True )) 
            if x.number_of_grid_points[2] > 1:
                if node_offset[2] == quantities.zero:
                    ghost_nodes.append(Node(code, mesh_length * (0, 0, 1) + node_offset * (1,1,0), True )) 
                if node_offset[2]  + node_size[2] == mesh_length[2]:
                    ghost_nodes.append(Node(code, node_size * (0, 0, -1) + node_offset * (1,1,0), True )) 
                raise Exception("3D not finished yet!")
                if node_offset[1] == quantities.zero:
                    ghost_nodes.append(Node(code, mesh_length * (0, 1, 0) + node_offset * (1,0,1), True )) 
                if node_offset[1] + node_size[1] == mesh_length[1]:
                    ghost_nodes.append(Node(code, node_size * (0, -1, 0)  + node_offset * (1,0,1), True ))  
                if node_offset[0] == quantities.zero and node_offset[1] == quantities.zero:
                    ghost_nodes.append(Node(code, mesh_length * (1, 1, 0) + node_offset * (0,0,1), True )) 
                if node_offset[0] == quantities.zero and node_offset[1] + node_size[1] == mesh_length[1]:
                    ghost_nodes.append(Node(code, mesh_length * (1, 0, 0) + (node_size * (0, -1, 0)) + node_offset * (0,0,1), True )) 
                if node_offset[1] == quantities.zero and node_offset[0] + node_size[0] == mesh_length[0]:
                    ghost_nodes.append(Node(code, mesh_length * (0, 1, 0) + (node_size * (-1, 0, 0)) + node_offset * (0,0,1), True )) 
                if node_offset[1] + node_size[1] == mesh_length[1] and node_offset[0] + node_size[0] == mesh_length[0]:
                    ghost_nodes.append(Node(code, node_size * (-1, -1, 0) + node_offset * (0,0,1), True )) 
        return ghost_nodes
        
    def stop(self):
        for x in self.nodes:
            if not x.is_live:
                x.stop()
                
    def get_hydro_state_at_point(self, x, y, z):
        rho_result = None
        rhovx_result = None
        rhovy_result = None
        rhovz_result = None
        rhoe_result = None
        for node in self.nodes:
            if not node.is_ghost:
                rho, rhovx, rhovy, rhovz, rhoe = node.code.get_hydro_state_at_point(
                    x - node.offset[0],
                    y - node.offset[1],
                    z - node.offset[2]
                )
                if rho_result is None:
                    rho_result = rho
                    rhovx_result = rhovx
                    rhovy_result = rhovy
                    rhovz_result = rhovz
                    rhoe_result = rhoe
                else:
                    rho_result += rho
                    rhovx_result += rhovx
                    rhovy_result += rhovy
                    rhovz_result += rhovz
                    rhoe_result += rhoe
        return (
            rho_result,
            rhovx_result,
            rhovy_result,
            rhovz_result,
            rhoe_result
        )
    
    def get_timestep(self):
        timesteps = [] | t_unit
        for x in self.nodes:
            if not x.is_ghost:
                timesteps.append(x.code.get_timestep())
        return timesteps.min()
            
        
