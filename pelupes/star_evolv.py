from amuse.support.units import units

from amuse.legacy.sse.interface import SSE
from amuse.legacy.evtwin.interface import EVtwin
from amuse.legacy.mesa.interface import MESA

from amuse.support.data import core

se=MESA()

star=core.Particle()
star.mass= 0.831 | units.MSun

se.initialize_module_with_current_parameters()

star=se.particles.add_particle(star)
se.initialize_stars()

stopped_evolving = False
while star.stellar_type.value_in(units.stellar_type) < 10 and not stopped_evolving:
    previous_age = star.age
    print star.age,star.stellar_type
    try:
        se.evolve_model()
        stopped_evolving = (star.age == previous_age) # Check whether the age has stopped increasing
    except Exception as ex:
        print str(ex)
        stopped_evolving = True
if stopped_evolving:
    print "Age did not increase during timestep. Aborted evolving..."
print " ... evolved model to t = " + str(star.age.as_quantity_in(units.Myr))
print "Star has now become a: ", star.stellar_type, "(stellar_type: "+str(star.stellar_type.value_in(units.stellar_type))+")"
print

