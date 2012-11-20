from amuse.community.evtwin.interface import EVtwin

from amuse.units import units

from amuse.datamodel import Particles

p=Particles(2,mass=[1,2]|units.MSun)

se=EVtwin()

se.particles.add_particles(p)

q=p[0].as_particle_in_set(se.particles)

q.evolve_one_step()
