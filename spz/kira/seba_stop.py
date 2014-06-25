from amuse.lab import *
from amuse.units.optparse import OptionParser
    
def main(M=1.0|units.MSun, z=0.02, model_time=4700|units.Myr):

    stellar = SeBa()#redirection="none")
    stellar.parameters.metallicity = z
    stars = Particles(3)
    stars[0].mass = 0.6*M
    stars[1].mass = M
    stars[2].mass = 0.5*M
    stellar.particles.add_particles(stars)
    stellar.commit_particles()
    stellar.stopping_conditions.supernova_detection.enable()

    L_init = stellar.particles.luminosity
    dt = 1 | units.Myr
    time = 0 | units.Myr
    while stellar.particles[0].age<model_time:
        stellar.evolve_model(time+dt)
        time = stellar.model_time
        print "time=", time, ", t=", stellar.particles.age, " type", stellar.particles.mass, stellar.particles.stellar_type.number
#        if stellar.stopping_conditions.supernova_detection.is_set():
#            print "Stellar supernova stopping condition is set"
#        print "star=", stellar.particles
        print "age=", stellar.particles.age
        
    import numpy 
    print "L=", numpy.log10(stellar.particles[0].luminosity.value_in(units.LSun))

    stellar.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-M", unit= units.MSun,
                      dest="M", type="float",default = 100.0 | units.MSun,
                      help="stellar mass [100.0] %unit")
    result.add_option("-t", unit = units.Myr,
                      dest="model_time", type="float", 
                      default = 4.5|units.Myr,
                      help="end time of the simulation [4.7] %unit")
    result.add_option("-z", dest="z", type="float", default = 0.02,
                      help="metalicity [0.02]")
    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", 
                      preferred_units = [units.MSun, units.parsec, units.Myr], 
                      precision = 16, prefix = "", 
                      separator = " [", suffix = "]")
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
