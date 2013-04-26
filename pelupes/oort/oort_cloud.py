"""
  oort cloud simulation
"""
import numpy
from amuse.lab import *
from amuse.couple import bridge
from amuse.community.gadget2.interface import Gadget2
from amuse.datamodel import ParticlesSuperset
from matplotlib import pyplot

from amuse.community.kepler2.interface import Kepler2

if __name__ in ('__main__', '__plot__'):
    t_end = 10**9 | units.yr

##### construct a system (an object in amuse/python)
# convert N-body system to Oort-cloud system
    converter= nbody_system.nbody_to_si(1.0|units.MSun, 10000. |units.AU)
    grav=Kepler2(converter)
    
##### creat particles for oort object
# comets
    comets = read_set_from_file('generate_id_6d.txt', 'txt', attribute_names=['x', 'y', 'z','vx', 'vy', 'vz' ],attribute_types=[units.AU, units.AU,  units.AU, units.AU/units.yr, units.AU/units.yr, units.AU/units.yr])
    N= len(comets.x)
    masses = numpy.zeros(N) | units.MSun
    comets.mass = masses
# sun
    sun = Particles(1)
    sun.mass = 1.0 | units.MSun
    sun.x = 0. | units.AU
    sun.y = 0. | units.AU
    sun.z = 0. | units.AU
    sun.vx = 0. | units.AU/units.yr
    sun.vy = 0. | units.AU/units.yr
    sun.vz = 0. | units.AU/units.yr
# combine all particles
    sun.add_particles(comets)

    grav.particles.add_particles(sun)

    rmin=(1-grav.orbiters_astro_centric.eccentricity)*grav.orbiters_astro_centric.semi_major_axis
    
    radius=15. | units.AU
    sel=numpy.where( rmin < radius)[0]
    
    for i in sel:
      print grav.orbiters_astro_centric[i].get_next_radial_crossing_time(radius).in_(units.Myr)
 
    model_time = 0 | units.yr

    times = 10**3*units.yr(range(0,99,1))
    f=pyplot.figure(figsize=(10,10))

    for i,ttarget in enumerate(times):
        grav.evolve_model(ttarget)
        print "time=", ttarget.in_(units.yr)
        x=grav.particles.x.value_in(units.AU)
        y=grav.particles.y.value_in(units.AU)
        """
        pyplot.plot(x,y,'r .')
        pyplot.plot([0.],[0.],'b +')
        pyplot.xlim(-40000,40000)
        pyplot.ylim(-40000,40000)
        """
    pyplot.show()        
    grav.stop()
