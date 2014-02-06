import os
import os.path
import shutil

from amuse.units import units
from amuse.datamodel import Particle

from amuse.ext.star_to_sph import pickle_stellar_model
#~from amuse.community.evtwin.interface import EVtwin as stellar_evolution_code
from amuse.community.mesa.interface import MESA as stellar_evolution_code


def to_model_directory():
    model_directory = os.path.join(os.getcwd(), "giant_models_"+stellar_evolution_code.__name__)
    if not os.path.exists(model_directory):
        os.mkdir(model_directory)
    shutil.copy(__file__, model_directory)
    os.chdir(model_directory)

def evolve_giant(giant, start_radius, stop_radius):
    stellar_evolution = stellar_evolution_code(redirection="file", redirect_file=stellar_evolution_code.__name__+".log")
    giant_in_code = stellar_evolution.particles.add_particle(giant)
    print "\nEvolving with", stellar_evolution_code.__name__
    while (giant_in_code.radius < start_radius):
        giant_in_code.evolve_one_step()
        print giant_in_code.radius, giant_in_code.age
    
    print "Stellar radius exceeds {0}, now saving model every step...".format(start_radius)
    print giant_in_code.as_set()
    
    i = 0
    while (giant_in_code.radius < stop_radius):
        giant_in_code.evolve_one_step()
        print giant_in_code.radius, giant_in_code.age
        pickle_file_name = "model_{0:=04}_{1:.1f}.pkl".format(i, giant_in_code.radius.value_in(units.RSun))
        pickle_stellar_model(giant_in_code, pickle_file_name)
        i += 1
    print "Finished: stellar radius {0} exceeds {1}, ".format(giant_in_code.radius, stop_radius)
    stellar_evolution.stop()


if __name__ == "__main__":
    to_model_directory()
    giant = Particle(mass = 10.0 | units.MSun)
    evolve_giant(giant, 800|units.RSun, 850|units.RSun)
