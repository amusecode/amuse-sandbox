import os
import sys
import numpy
from matplotlib import pyplot
import matplotlib.cm as cm


from amuse.support.data import core
from amuse.ext import cloud
from amuse.support.units import generic_unit_system
from amuse.community.athena.interface import Athena

def new_cloudshock_code(N = 32):
    L = 10.0
    grid = core.Grid.create((N,4*N,N), [10.0, 40.0, 10.0] | generic_unit_system.length)
    core.Grid.add_global_vector_attribute("position", ["x","y","z"])
    
    cloud.fill_grid_with_cloud_shock(
        grid, 
        center = [5.0, 5.0, 5.0] | generic_unit_system.length,
        radius = 1.0 | generic_unit_system.length,
    )
    
    
    instance=Athena()
    instance.initialize_code()
    instance.parameters.gamma = 5.0 / 3.0
    instance.parameters.courant_number = 0.3
    instance.set_boundary("reflective","reflective","outflow","outflow",
                          "reflective","reflective")
    instance.setup_mesh(N, 4*N, N, L,4*L,L)
    instance.commit_parameters()
    
    channel = grid.new_channel_to(instance.grid)
    channel.copy()
    instance.initialize_grid()
    return instance

if __name__=="__main__":
    N=32
    instance=new_cloudshock_code(N)
    
    rc=1.
    gamma=5./3
    xi=10.
    cs=numpy.sqrt((gamma-1))
    cs_out=numpy.sqrt((gamma-1)*xi)
    vs=cs_out*2.7
    tau=1.6*2*rc*xi**0.5/vs
    
    print tau
    
    pyplot.figure(figsize=(12,12))
    
    instance.evolve(0.00*tau | generic_unit_system.time ) 
    print instance.get_time()
    instance.evolve(0.5*tau | generic_unit_system.time ) 
    print instance.get_time()
    
    rho =instance.grid.rho[...,...,N/2]
    pyplot.subplot(1,1,1)
    print "bla"
    print rho[0]
    pyplot.imshow( rho.value_in(generic_unit_system.mass/generic_unit_system.length**3))
    
    pyplot.savefig("cloud.png")
