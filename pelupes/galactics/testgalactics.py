from amuse.community.galactics.interface import GalactICs, GalactICsInterface

from amuse.units import units,nbody_system

from matplotlib import pyplot

convert = nbody_system.nbody_to_si(1.0 | units.kpc, 1.0e6 | units.MSun)

instance = GalactICs(unit_converter=convert)

instance.parameters.halo_random_seed=12345
instance.parameters.bulge_number_of_particles=1000
instance.parameters.halo_number_of_particles=5000
instance.parameters.disk_number_of_particles=1000
instance.parameters.output_directory="./"
print instance.parameters

print instance.model_present()

instance.commit_parameters()

instance.generate_particles()

p=instance.particles
halo=instance.halo_particles
bulge=instance.bulge_particles
disk=instance.disk_particles

print halo

print len(p)

print p.total_mass().in_(units.MSun)

pyplot.plot(disk.x.value_in(units.kpc),disk.z.value_in(units.kpc),'r.')
pyplot.xlim(-30,30)
pyplot.ylim(-30,30)
pyplot.show()
