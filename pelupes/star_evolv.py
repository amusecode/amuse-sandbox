from amuse.community.sse.interface import SSE
from amuse.community.evtwin.interface import EVtwin
from amuse.community.mesa.interface import MESA

from amuse.units import units

from amuse import datamodel
se=SSE()
se.initialize_code()

star=datamodel.Particle()
star.mass= 0.831 | units.MSun

se.commit_parameters()

star=se.particles.add_particle(star)
se.commit_particles()

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

