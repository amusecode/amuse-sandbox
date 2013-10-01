import numpy
import matplotlib
matplotlib.use("Agg")

from amuse.units import units, nbody_system, constants
#~from amuse.units.quantities import zero
from amuse.datamodel import Particles, Particle
#~from amuse.support.exceptions import AmuseException
from amuse.plot import pynbody_column_density_plot, HAS_PYNBODY

#~from amuse.community.bhtree.interface import BHTree
from amuse.community.hermite0.interface import Hermite
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
#~from amuse.community.seba.interface import SeBa
#~from amuse.community.evtwin.interface import EVtwin
from amuse.community.mesa.interface import MESA

from amuse.couple.collision_handler import CollisionHandler
from amuse.ext.plotting_hydro import new_plotting_hydrodynamics_code
from amuse.ext.hydro_collision import StellarEncounterInHydrodynamics


def new_colliders(masses, separation):
    colliders = Particles(2)
    colliders.mass = masses * 1.0
    colliders.position = [0.0, 0.0, 0.0] | units.RSun
    colliders.velocity = [0.0, 0.0, 0.0] | units.km / units.s
    colliders[0].x += separation
    colliders[0].vy += 2 * (constants.G * colliders.total_mass() / separation).sqrt()
    colliders.move_to_center()
    return colliders

def new_gravity(particles):
    convert_nbody = nbody_system.nbody_to_si(particles.total_mass(), 10.0 | units.RSun)
    gravity = Hermite(convert_nbody)
    gravity.particles.add_particles(particles)
    return gravity

def new_stellar(particles):
    stellar = MESA(redirection="file", redirect_file="mesa_out.log")
    stellar.particles.add_particles(particles)
    stellar.evolve_model(3.0 | units.Myr)
    return stellar

def new_collision_handler(number_of_sph_particles, stellar, gravity):
    max_density = max(stellar.particles[0].get_density_profile().amax(), stellar.particles[1].get_density_profile().amax())
    density_resolution = 0.1 * stellar.particles.mass.sum() / (number_of_sph_particles | units.RSun**2)
    density_min = numpy.log10(density_resolution.value_in(constants.proton_mass.to_unit() / units.cm**2))
    density_max = numpy.log10((max_density * (10|units.RSun)).value_in(constants.proton_mass.to_unit() / units.cm**2))
    
    dt = 0.01 * ((1.0|units.RSun**3 / units.MSun) / constants.G).sqrt()
    collision_code = StellarEncounterInHydrodynamics(
        number_of_sph_particles,
        new_plotting_hydrodynamics_code(
#~            Gadget2, 
            Fi, 
            0.02|units.hour, 
            plot_function = pynbody_column_density_plot,
            plot_function_arguments = dict(width=200|units.RSun, vmin=density_min, vmax=density_max)
        ), 
        initial_separation = 2,
#~        hydrodynamics_arguments = dict(redirection="file", redirect_file="gadget.out"),
        hydrodynamics_arguments = dict(redirection="file", redirect_file="fi.out"),
#~        hydrodynamics_parameters = dict(time_limit_cpu=1|units.day),
        hydrodynamics_parameters = dict(timestep=dt, eps_is_h_flag=True),
        verbose = True,
        debug = False
    )
    collision_handler = CollisionHandler(
        collision_code, 
        stellar_evolution_code=stellar, 
        gravity_code=gravity
    )
    return collision_handler


def run(number_of_sph_particles = 100000, masses = [40, 20]|units.MSun, separation = 0.5|units.RSun):
    print "Simulating collision with {0} SPH particles".format(number_of_sph_particles)
    colliders = new_colliders(masses, separation)
    gravity = new_gravity(colliders)
    stellar = new_stellar(colliders)
    collision_handler = new_collision_handler(number_of_sph_particles, stellar, gravity)
    
    result = collision_handler.handle_collision(colliders[0], colliders[1])
    
    gravity.stop()
    stellar.evolve_model(0.02 | units.Myr)
    print stellar.particles
    stellar.stop()
    return result.mass


if __name__ == "__main__":
    number_of_sph_particles = [5000, 10000, 20000, 50000, 100000, 200000, 500000]

#~    number_of_sph_particles = [10000, 20000]
#~    [0.5, 1, 3, 5]
#~[10000, 20000] [[quantity<[35.016179016] 1.98892e+30 * kg>, quantity<[20.6873650802, 11.8748891614] 1.98892e+30 * kg>, quantity<[36.1781614105, 14.227400147] 1.98892e+30 * kg>, quantity<[39.2199405654, 19.3939086747] 1.98892e+30 * kg>], [quantity<[35.5467636472] 1.98892e+30 * kg>, quantity<[32.6307745706] 1.98892e+30 * kg>, quantity<[36.0337547055, 14.3512639382] 1.98892e+30 * kg>, quantity<[39.2888635366, 19.4959986343] 1.98892e+30 * kg>]]
    
    mass = []
    for n in number_of_sph_particles:
        new_list = []
#~        for sep in [0.5, 1, 3, 5] | units.RSun: # How to do head-on (vy should be ~0, but I have to think about vx)?
        for sep in [1.0] | units.RSun: # How to do head-on (vy should be ~0, but I have to think about vx)?
            try:
                result = run(number_of_sph_particles=n, separation=sep)
            except Exception as result:
                pass
            new_list.append(result)
        mass.append(new_list)
    
    print number_of_sph_particles, mass
