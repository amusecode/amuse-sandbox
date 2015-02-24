import multiplexing_gd

from amuse.test.amusetest import TestWithMPI

if __name__ == "__main__":
    x = multiplexing_gd.MultiplexingGravitationalDynamicsCode()
    print x.particles
from amuse import datamodel


from amuse.units import nbody_system


from amuse.units import units


import numpy


class TestGravitationalDynamicsInterface(TestWithMPI):
    
    def test01(self):
        x = multiplexing_gd.MultiplexingGravitationalDynamicsInterface()
        index, error = x.new_multiplexed_particle(1, 3, 4, 5, 6, 7, 8, 9, 1)
        self.assertEquals(error, 0)
        self.assertEquals(index, 1)
        mass, error = x.get_mass(index)
        self.assertEquals(error, 0)
        self.assertEquals(mass, 1)
        x, y, z, error = x.get_position(index)
        self.assertEquals(error, 0)
        self.assertEquals(x, 3)
        self.assertEquals(y, 4)
        self.assertEquals(z, 5)
    def test02(self):
        x = multiplexing_gd.MultiplexingGravitationalDynamicsInterface()
        index, error = x.new_multiplexed_particle([1,12], [3,13], [4,14], [5,15], [6,16], [7,17], [8,18], [9,19], [1,1])
        self.assertEquals(error, 0)
        self.assertEquals(index, [1,2])
        mass, error = x.get_mass(index)
        self.assertEquals(error, 0)
        self.assertEquals(mass, [1,12])
        x, y, z, error = x.get_position(index)
        self.assertEquals(error, 0)
        self.assertEquals(x, [3,13])
        self.assertEquals(y, [4,14])
        self.assertEquals(z, [5,15])

    def test03(self):
        x = multiplexing_gd.MultiplexingGravitationalDynamicsInterface()
        index, error = x.new_multiplexed_particle([1,12], [3,13], [4,14], [5,15], [6,16], [7,17], [8,18], [9,19], [1,1])
        self.assertEquals(error, 0)
        self.assertEquals(index, [1,2])
        error = x.set_mass(index, [21,212])
        print "error:", error
        self.assertEquals(error, 0)
        mass, error = x.get_mass(index)
        self.assertEquals(error, 0)
        self.assertEquals(mass, [21,212])
        error = x.set_position(index, [23,213], [24, 214], [25, 215])
        self.assertEquals(error, 0)
        x, y, z, error = x.get_position(index)
        self.assertEquals(error, 0)
        self.assertEquals(x, [23,213])
        self.assertEquals(y, [24,214])
        self.assertEquals(z, [25,215])
    def test04(self):
        x = multiplexing_gd.MultiplexingGravitationalDynamicsInterface()
        index, error = x.new_multiplexed_particle([1,12], [3,13], [4,14], [5,15], [6,16], [7,17], [8,18], [9,19], [1,1])
        self.assertEquals(error, 0)
        self.assertEquals(index, [1,2])
        error = x.delete_particle(index[0])
        self.assertEquals(error, 0)
        


class TestMultiplexingGravitationalDynamics(TestWithMPI):
    
    def test01(self):
        x = multiplexing_gd.MultiplexingGravitationalDynamicsCode()
        particles = datamodel.Particles(5)
        particles.mass = 1.0 | nbody_system.mass
        particles.position = [0,0,0] | nbody_system.length
        particles.velocity = [0,0,0] | nbody_system.speed
        #particles.index_of_the_set = 0
        print particles
        x.particles.add_particles(particles)
        print x.particles
        self.assertEquals(len(x.particles), 5)
        self.assertEquals(x.particles[0].mass, 1.0 | nbody_system.mass)
    def test02(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        hermite = multiplexing_gd.MultiplexingGravitationalDynamicsCode(convert_nbody)
        
        try:
            hermite.parameters.epsilon_squared = 0.0 | units.AU**2
            #hermite.parameters.end_time_accuracy_factor = 0.0

            stars = self.new_system_of_sun_and_earth()
            earth = stars[1]

            hermite.particles.add_particles(stars)

            hermite.evolve_model(365.0 | units.day)
            hermite.particles.copy_values_of_all_attributes_to(stars)

            position_at_start = earth.position.value_in(units.AU)[0]
            position_after_full_rotation = earth.position.value_in(units.AU)[0]
            self.assertAlmostRelativeEqual(position_at_start, position_after_full_rotation, 6)

            hermite.evolve_model(365.0 + (365.0 / 2) | units.day)

            hermite.particles.copy_values_of_all_attributes_to(stars)
            position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
            self.assertAlmostRelativeEqual(-position_at_start, position_after_half_a_rotation, 3)

            hermite.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)

            hermite.particles.copy_values_of_all_attributes_to(stars)
            position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
            self.assertAlmostRelativeEqual(-position_at_start, position_after_half_a_rotation, 3)

            hermite.cleanup_code()
        finally:
            hermite.stop()
    def new_system_of_sun_and_earth(self):
        stars = datamodel.Stars(2)
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = units.m(numpy.array((0.0,0.0,0.0)))
        sun.velocity = units.ms(numpy.array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)

        earth = stars[1]
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371)     
        earth.position = units.km(numpy.array((149.5e6,0.0,0.0)))
        earth.velocity = units.ms(numpy.array((0.0,29800,0.0)))
        
        return stars
        

    def test03(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        hermite = multiplexing_gd.MultiplexingGravitationalDynamicsCode(convert_nbody)
        
        try:
            hermite.parameters.epsilon_squared = 0.0 | units.AU**2
            #hermite.parameters.end_time_accuracy_factor = 0.0

            stars = self.new_system_of_sun_and_earth()
            stars2 = self.new_system_of_sun_and_earth()
            earth = stars[1]

            hermite.particles.add_particles(stars)
            stars2.index_of_the_set = 1
            hermite.particles.add_particles(stars2)
            
            hermite.evolve_model(365.0 | units.day)
            hermite.particles.copy_values_of_all_attributes_to(stars)

            position_at_start = earth.position.value_in(units.AU)[0]
            position_after_full_rotation = earth.position.value_in(units.AU)[0]
            self.assertAlmostRelativeEqual(position_at_start, position_after_full_rotation, 6)

            hermite.evolve_model(365.0 + (365.0 / 2) | units.day)

            hermite.particles.copy_values_of_all_attributes_to(stars)
            position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
            self.assertAlmostRelativeEqual(-position_at_start, position_after_half_a_rotation, 3)

            hermite.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)

            hermite.particles.copy_values_of_all_attributes_to(stars)
            position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
            self.assertAlmostRelativeEqual(-position_at_start, position_after_half_a_rotation, 3)

            hermite.cleanup_code()
        finally:
            hermite.stop()

    def test04(self):
        #convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        hermite = multiplexing_gd.MultiplexingMercuryCode()
        hermite.set_code("Mercury")
        try:
            hermite.parameters.epsilon_squared = 0.0 | units.AU**2
            #hermite.parameters.end_time_accuracy_factor = 0.0

            stars = self.new_system_of_sun_and_earth()
            earth = stars[1]

            hermite.particles.add_particles(stars)
            stars2 = self.new_system_of_sun_and_earth()
            stars2[1].mass *= 2.0
            stars2.index_of_the_set = 1
            hermite.particles.add_particles(stars2)
            
            hermite.evolve_model(365.0 | units.day)
            hermite.particles.copy_values_of_all_attributes_to(stars)
            
            position_at_start = earth.position.value_in(units.AU)[0]
            position_after_full_rotation = earth.position.value_in(units.AU)[0]
            self.assertAlmostRelativeEqual(position_at_start, position_after_full_rotation, 6)

            hermite.evolve_model(365.0 + (365.0 / 2) | units.day)

            hermite.particles.copy_values_of_all_attributes_to(stars)
            position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
            self.assertAlmostRelativeEqual(-position_at_start, position_after_half_a_rotation, 2)

            hermite.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)

            hermite.particles.copy_values_of_all_attributes_to(stars)
            position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
            self.assertAlmostRelativeEqual(-position_at_start, position_after_half_a_rotation, 2)

            hermite.cleanup_code()
        finally:
            hermite.stop()





