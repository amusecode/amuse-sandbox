"""

"""


import numpy
from matplotlib import pyplot
from optparse import OptionParser


from amuse.community.fi.interface import Fi
from amuse.community.mocassin.interface import Mocassin


from amuse.units import units
from amuse.units import nbody_system

from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.datamodel import Particles
from amuse.datamodel import Grid


def new_grid_from_code(sph, number_of_cells=100, box_size=1 | units.AU):
    grid = Grid.create(
        (number_of_cells,number_of_cells,number_of_cells),
        ([1,1,1]|units.none) * box_size
    )
    
    half_length = box_size / 2 
    grid.position -= ([1,1,1]|units.none) * half_length
    
    grid.vx = 0 | units.kms
    grid.vy = 0 | units.kms
    grid.vz = 0 | units.kms
    
    rho,_,_,_,_ = sph.get_hydro_state_at_point(grid.x.flatten(),grid.y.flatten(),grid.z.flatten(),grid.vx.flatten(),grid.vy.flatten(),grid.vz.flatten())
    
    grid.rho = rho.reshape(grid.shape)
    
    return grid
    

def main(number_of_gas_particles = 2000, t_end = 1, star_mass = 1):
    #numpy.random.seed(1234)
    
    t_end = t_end | units.yr
    star_mass = star_mass | units.MSun
        
    convert=nbody_system.nbody_to_si(star_mass, 1. | units.AU)
    proto=ProtoPlanetaryDisk(number_of_gas_particles,convert_nbody=convert,densitypower=1.5,Rmin=4,Rmax=20,q_out=1.)
    gas=proto.result
    gas.h_smooth=0.06 | units.AU
         
    sun=Particles(1)
    sun.mass=star_mass
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
        
    #sph.evolve_model(t_end)    
    number_of_cells = 41
    box_size = 50 | units.AU
    print 1
    grid=new_grid_from_code(sph,number_of_cells=number_of_cells,box_size = box_size)
    print 2
    
    sph.stop()
    
    radiative_transfer = Mocassin(redirection = "none", number_of_workers = 3)
    radiative_transfer.set_input_directory(radiative_transfer.get_default_input_directory())
    radiative_transfer.initialize_code()
    radiative_transfer.set_symmetricXYZ(False)    
    radiative_transfer.parameters.length_x = box_size
    radiative_transfer.parameters.length_y = box_size
    radiative_transfer.parameters.length_z = box_size
    radiative_transfer.parameters.mesh_size = (number_of_cells, number_of_cells, number_of_cells)
    radiative_transfer.commit_parameters()
    
    sun.temperature = 20000 | units.K
    sun.luminosity = 1.0  | units.LSun
    radiative_transfer.grid.hydrogen_density = (grid.rho / (1.0 | units.amu))
    
    radiative_transfer.commit_grid()
    radiative_transfer.particles.add_particle(sun)
    radiative_transfer.commit_particles()
    
    radiative_transfer.step()
    
    etemp = radiative_transfer.grid.electron_temperature[...,...,number_of_cells / 2 + 1]
    print etemp[..., number_of_cells / 2 + 1]
    pyplot.figure(figsize=(8,8))
    pyplot.imshow(
        etemp.value_in(units.K),
        extent = (([-0.5,0.5,-0.5,0.5]|units.none) * box_size).value_in(units.AU),
    )
    pyplot.title(t_end)
    pyplot.xlabel('AU')
    pyplot.ylabel('AU')
    pyplot.show()
         
def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-p", "--ngasparticles", 
        default = 20000,
        dest="number_of_gas_particles",
        help="number of gas particles in the disk",
        type="int"
    )
    result.add_option(
        "-t", "--endtime", 
        default = 1,
        dest="t_end",
        help="number of years to evolve the protoplanetary disk",
        type="float"
    )
    result.add_option(
        "-m", "--starmass", 
        default = 1,
        dest="star_mass",
        help="mass of the central star, in solar masses",
        type="float"
    )
    return result
         
if __name__ == "__plot__":
    main(20000, 1, 1)
    
if __name__ == "__main__":
    print units.amu.value_in(units.kg)
    options, arguments = new_option_parser().parse_args()
    main(**options.__dict__)

