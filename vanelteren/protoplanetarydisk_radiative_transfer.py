"""

"""


import numpy

from matplotlib import pyplot

from amuse.community.fi.interface import Fi
from amuse.community.mocassin.interface import Mocassin


from amuse.units import units
from amuse.units import nbody_system

from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.datamodel import Particles
from amuse.datamodel import Grid


def make_grid(sph,N=100,L=1):
    grid = Grid.create((N,N,N), [L, L, L] | units.AU)
    
    half_length = L / 2 
    grid.position -= (half_length,half_length,half_length)  | units.AU
    
    grid.vx = 0 | units.kms
    grid.vy = 0 | units.kms
    grid.vz = 0 | units.kms
    
    rho,_,_,_,_ = sph.get_hydro_state_at_point(grid.x.flatten(),grid.y.flatten(),grid.z.flatten(),grid.vx.flatten(),grid.vy.flatten(),grid.vz.flatten())
    
    grid.rho = rho.reshape(grid.shape)
    
    return grid
    
if __name__ in ("__main__","__plot__"):

    N=20000
    tend=1. | units.yr
    Mstar=1. | units.MSun
        
    convert=nbody_system.nbody_to_si(Mstar, 1. | units.AU)
    proto=ProtoPlanetaryDisk(N,convert_nbody=convert,densitypower=1.5,Rmin=4,Rmax=20,q_out=1.)
    gas=proto.result
    gas.h_smooth=0.06 | units.AU
         
    sun=Particles(1)
    sun.mass=Mstar
    sun.radius=2. | units.AU
    sun.x=0.|units.AU
    sun.y=0.|units.AU
    sun.z=0.|units.AU
    sun.vx=0.|units.kms
    sun.vy=0.|units.kms
    sun.vz=0.|units.kms
 
    sph=Fi(convert)
 
    sph.parameters.use_hydro_flag=True
    sph.parameters.radiation_flag=False
    sph.parameters.self_gravity_flag=True
    sph.parameters.gamma=1.
    sph.parameters.isothermal_flag=True
    sph.parameters.integrate_entropy_flag=False
    sph.parameters.timestep=0.125 | units.yr  

    sph.gas_particles.add_particles(gas)
    sph.particles.add_particles(sun)
        
    #sph.evolve_model(tend)    
            
    L=50
    grid_size = 41
    print 1
    grid=make_grid(sph,N=grid_size,L=L)
    print 2
    
    sph.stop()
    
    radiative_transfer = Mocassin() #redirection = "none")
    radiative_transfer.set_input_directory(radiative_transfer.get_default_input_directory())
    radiative_transfer.initialize_code()
    radiative_transfer.set_symmetricXYZ(False)
    radiative_transfer.set_constant_hydrogen_density(100.0 | (1/units.cm**3))
    
    radiative_transfer.parameters.length_x = L | units.AU
    radiative_transfer.parameters.length_y = L | units.AU
    radiative_transfer.parameters.length_z = L | units.AU
    radiative_transfer.parameters.mesh_size = (grid_size, grid_size, grid_size)
    
    radiative_transfer.commit_parameters()
    sun.temperature = 20000 | units.K
    sun.luminosity = 1.0  | units.LSun
    print radiative_transfer.grid.shape
    print radiative_transfer.grid.z[21]
    sys.exit(0)
    print 3
    #radiative_transfer.set_has_constant_hydrogen_density(True)
    radiative_transfer.commit_grid()
    radiative_transfer.particles.add_particle(sun)
    radiative_transfer.commit_particles()
    
    rho = (grid.rho * (1.0 / (1.0 | units.kg)))
    print radiative_transfer.grid.shape
    radiative_transfer.grid.hydrogen_density = rho
    print radiative_transfer.grid.x[5]
    print "done"
    hden = radiative_transfer.grid.hydrogen_density[...,...,20]
    print     hden
    print rho[...,...,5]
    pyplot.figure(figsize=(8,8))
    pyplot.imshow(hden.value_in(1/units.cm**3),
        extent=[-L/2,L/2,-L/2,L/2])#,vmin=10,vmax=15)    
    pyplot.title(tend)
    pyplot.xlabel('AU')
    pyplot.show()
         
