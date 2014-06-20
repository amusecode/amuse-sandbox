from amuse.community.galactics.gas_interface import GaslactICs, GaslactICsInterface

from amuse.units import units,nbody_system

from matplotlib import pyplot

convert = nbody_system.nbody_to_si(1.0 | units.kpc, 1.0e9 | units.MSun)

print convert.to_nbody(8 | units.kms)

instance = GaslactICs(unit_converter=convert)

#instance.parameters.halo_random_seed=12345
#instance.parameters.bulge_number_of_particles=1000
#instance.parameters.halo_number_of_particles=50
#instance.parameters.disk_number_of_particles=1000
instance.parameters.output_directory="./"
#print instance.parameters

instance.commit_parameters()

instance.generate_particles()

p=instance.particles
halo=instance.halo_particles
bulge=instance.bulge_particles
disk=instance.disk_particles
gas=instance.gas_particles

print len(p)

print p.total_mass().in_(units.MSun)

print gas.u[0]**0.5

#pyplot.plot(halo.x.value_in(units.kpc),halo.z.value_in(units.kpc),'y.')
#pyplot.plot(bulge.x.value_in(units.kpc),bulge.z.value_in(units.kpc),'b.')
#pyplot.plot(disk.x.value_in(units.kpc),disk.z.value_in(units.kpc),'g.')
#pyplot.plot(gas.x.value_in(units.kpc),gas.z.value_in(units.kpc),'r.')
#pyplot.xlim(-30,30)
#pyplot.ylim(-30,30)

r=(gas.x**2+gas.y**2)**0.5
pyplot.plot(r.value_in(units.kpc),gas.z.value_in(units.kpc),'r.')
pyplot.xlim(0,30)
pyplot.ylim(-5,5)

pyplot.show()
