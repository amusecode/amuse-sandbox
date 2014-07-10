"""
    simple coupling of a stellar evolution and a gravitational code to simulate a cluster of stars.
"""

from amuse.datamodel.particles import AbstractParticleSet
from amuse.units.optparse import OptionParser
from amuse.units import units, quantities, nbody_system

from amuse.community.hermite0.interface import Hermite
from amuse.community.seba.interface import SeBa

from amuse.ic.plummer import new_plummer_model
from amuse.ic.salpeter import new_salpeter_mass_distribution

from matplotlib import pyplot
from amuse import plot as aplot

from amuse.support.console import set_printing_strategy

set_printing_strategy("custom", preferred_units = [units.MSun, units.parsec, units.Myr], precision=3)

class MultiParticlesOverlay(AbstractParticleSet):
    """An overlay of one or more particles set(s). The overlay
    provides a single interface to the particles in the multiple sets
    and stores extra attributes in the overlayed set

    >>> p1 = Particles(3)
    >>> p1.mass = [10.0, 20.0, 30.0] | units.kg
    >>> p2 = ParticlesOverlay(p1)
    >>> p2.radius = [4.0, 5.0, 6.0] | units.m
    >>> print len(p2)
    3
    >>> print p2.mass
    [10.0, 20.0, 30.0] kg
    """

    def __init__(self, particles, overlay_set = None):
        AbstractParticleSet.__init__(self)

        self._private.base_set = particles.as_set()
        if overlay_set is None:
            overlay_set = Particles(keys = self._private.base_set.key)
        self._private.overlay_set = overlay_set


    def can_extend_attributes(self):
        return self._private.overlay_set.can_extend_attributes()

    def __len__(self):
        return len(self._private.overlay_set)

    def _get_subsets_version(self):
        versions = [[x._get_version()] for x in self._private.particle_sets]
        return numpy.sum(versions)

    def _get_version(self):
        return self._private.overlay_set._get_version() +  self._private.base_set._get_version()


    def __getitem__(self, index):
        keys = self.get_all_keys_in_store()[index]

        if hasattr(keys, '__iter__'):
            return self._subset(keys)
        else:
            return Particle(keys, self)

    def _split_attributes(self, attributes):
        inbase = set(self._private.base_set.get_attribute_names_defined_in_store())
        attributes_inbase = []
        attributes_inoverlay = []
        indices_inbase = []
        indices_inoverlay = []
        for i,x in enumerate(attributes):
            if x in inbase:
                attributes_inbase.append(x)
                indices_inbase.append(i)
            else:
                attributes_inoverlay.append(x)
                indices_inoverlay.append(i)
        return (attributes_inbase, indices_inbase), (attributes_inoverlay, indices_inoverlay)

    def _split_attributes_and_values(self, attributes, values):
        inbase = set(self._private.base_set.get_attribute_names_defined_in_store())
        attributes_inbase = []
        attributes_inoverlay = []
        values_inbase = []
        values_inoverlay = []
        for x,y in zip(attributes, values):
            if x in inbase:
                attributes_inbase.append(x)
                values_inbase.append(y)
            else:
                attributes_inoverlay.append(x)
                values_inoverlay.append(y)
        return (attributes_inbase, values_inbase), (attributes_inoverlay, values_inoverlay)

    def add_particles_to_store(self, keys, attributes = [], values = []):
        (
            (attributes_inbase, values_inbase),
            (attributes_inoverlay, values_inoverlay)
        ) = self._split_attributes_and_values(attributes, values)


        self._private.base_set.add_particles_to_store(keys, attributes_inbase, values_inbase)
        self._private.overlay_set.add_particles_to_store(keys, attributes_inoverlay, values_inoverlay)

    def remove_particles_from_store(self, indices):
        indices = numpy.asarray(indices)

        self._private.base_set.remove_particles_from_store(indices[...,0])
        self._private.overlay_set.remove_particles_from_store(indices[...,1])



    def get_values_in_store(self, indices, attributes):

        (
            (attributes_inbase, indices_inbase),
            (attributes_inoverlay, indices_inoverlay)
        ) = self._split_attributes(attributes)


        if indices is None:
            indices0 = self._private.base_set.get_all_indices_in_store()
            indices1 = self._private.overlay_set.get_all_indices_in_store()
        else:
            indices0 = []
            indices1 = []
            for i0, i1 in indices:
                indices0.append(i0)
                indices1.append(i1)
            indices0 = numpy.asarray(indices0, dtype='int64')
            indices1 = numpy.asarray(indices1, dtype='int64')

        result = [None] * len(attributes)
        if len(attributes_inbase) > 0:
            values_inbase = self._private.base_set.get_values_in_store(indices0, attributes_inbase)
            for i, value in zip(indices_inbase, values_inbase):
                result[i] = value

        if len(attributes_inoverlay) > 0:
            values_inoverlay = self._private.overlay_set.get_values_in_store(indices1, attributes_inoverlay)
            for i, value in zip(indices_inoverlay, values_inoverlay):
                result[i] = value

        return result

    def set_values_in_store(self, indices, attributes, values):
        (
            (attributes_inbase, values_inbase),
            (attributes_inoverlay, values_inoverlay)
        ) = self._split_attributes_and_values(attributes, values)


        if indices is None:
            indices0 = self._private.base_set.get_all_indices_in_store()
            indices1 = self._private.overlay_set.get_all_indices_in_store()
        else:
            indices0 = []
            indices1 = []
            for i0, i1 in indices:
                indices0.append(i0)
                indices1.append(i1)
            indices0 = numpy.asarray(indices0, dtype='int64')
            indices1 = numpy.asarray(indices1, dtype='int64')

        if len(attributes_inbase) > 0:
            self._private.base_set.set_values_in_store(indices0, attributes_inbase, values_inbase)
        if len(attributes_inoverlay) > 0:
            self._private.overlay_set.set_values_in_store(indices1, attributes_inoverlay, values_inoverlay)


    def get_attribute_names_defined_in_store(self):
        result = list(self._private.base_set.get_attribute_names_defined_in_store())
        result.extend(self._private.overlay_set.get_attribute_names_defined_in_store())
        return result

    def get_all_keys_in_store(self):
        return self._private.overlay_set.get_all_keys_in_store()

    def get_all_indices_in_store(self):
        indices0 = self._private.base_set.get_all_indices_in_store()
        indices1 = self._private.overlay_set.get_all_indices_in_store()

        return zip(indices0, indices1)

    def get_indices_of_keys(self, keys):
        indices0 = self._private.base_set.get_indices_of_keys(keys)
        indices1 = self._private.overlay_set.get_indices_of_keys(keys)

        return zip(indices0, indices1)

    def has_key_in_store(self, key):
        return self._private.overlay_set.has_key_in_store(key)

    def _original_set(self):
        return self

def create_stars(number_of_stars, size):
    masses = new_salpeter_mass_distribution(number_of_stars, mass_min = 2|units.MSun)
    converter = nbody_system.nbody_to_si(masses.sum(), size)
    stars = new_plummer_model(number_of_stars, convert_nbody=converter)
    stars.mass = masses
    stars.zams_mass = masses

    return stars, converter

def output_results(stars, time):
    mass_loss = stars.zams_mass - stars.mass

    aplot.plot(stars.x, stars.y, "*")

    for x, y, mass_loss in zip(stars.x.number, stars.y.number, mass_loss):
        pyplot.annotate("%0.2f"%abs(mass_loss.number), xy=(x-1e17,y+5e16))

    pyplot.xlim([-2e18, 2e18])
    pyplot.ylim([-2e18, 2e18])

    name = "plots/plot_{0:=05}.png".format(int(time.number))
    print "creating", name
    pyplot.savefig(name)
    pyplot.close()

def my_gravity_and_stellar_evolution(number_of_stars, size, end_time, sync_timestep=1|units.Myr, output_timestep=2|units.Myr):
    stars, converter = create_stars(number_of_stars, size)

    stars = MultiParticlesOverlay(particles=stars)

    gravity = Hermite(converter)
    stellar = SeBa()

    stars.include_particles(gravity.particles)
    stars.include_particles(stellar.particles, source_of=["mass", "radius"])

    while gravity.model_time <= end_time:
        stars.sync()

        time = gravity.model_time + sync_timestep
        gravity.evolve_model(time)
        stellar.evolve_model(time)

        if time % plot_timestep < sync_timestep:
            output_results(stars, time)

def gravity_and_stellar_evolution(number_of_stars, size, end_time, sync_timestep=1|units.Myr, output_timestep=2|units.Myr):

    stars, converter = create_stars(number_of_stars, size)

    gravity = Hermite(converter)
    stellar = SeBa()

    gravity.particles.add_particles(stars)
    stellar.particles.add_particles(stars)

    channels = []
    channels.append(stellar.particles.new_channel_to(grav.particles, attributes=["mass", "radius"]))
    channels.append(stellar.particles.new_channel_to(stars))
    channels.append(grav.particles.new_channel_to(stars))

    while gravity.model_time <= end_time:
        for channel in channels:
            channel.copy()

        gravity.evolve_model(time)
        stellar.evolve_model(time)

        if time % plot_timestep < sync_timestep:
            output_results(stars, time)

def parse_arguments():
    parser = OptionParser()
    parser.add_option("-N", dest="number_of_stars", type="int", default=100,
        help="The number of stars in the cluster [%default].")
    parser.add_option("-s", dest="size", unit=units.parsec, default=10,
        help="The total size of the cluster [%default %unit].")
    parser.add_option("-t", dest="end_time", unit=units.Gyr, default=1,
        help="The end time of the simulation [%default %unit].")

    options, args = parser.parse_args()
    return options.__dict__

if __name__ == "__main__":
    options = parse_arguments()
    gravity_and_stellar_evolution(**options)
