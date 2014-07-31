"""
simple example to use kepler_orbiters
  integrates the solar system planets motion around the Sun
"""

from amuse.units import units, nbody_system, quantities
from amuse.datamodel import Particles
from amuse.community.kepler_orbiters.interface import Kepler
from matplotlib import pyplot

def sun_and_planets():
  """
  solar system on 2014-Jul-09 (JD 2456847.5)
    Sun + 8 planets
    http://ssd.jpl.nasa.gov/horizons.cgi
  """
  sun = Particles(1)
  sun[0].mass = 1.0 | units.MSun
  sun[0].position = (2.037427757106262E-03,
                     -1.670788401563819E-03,
                     -1.177108759951928E-04) | units.AU
  sun[0].velocity = (5.088355581933507E-06,
                     3.932938906525153E-06,
                     -1.245782562840000E-07) | (units.AU/units.day)
  sun[0].radius = 1.0 | units.RSun
  
  planets = Particles(8)
  mercury = planets[0]
  mercury.mass = 3.302e23 | units.kg
  mercury.position = (3.387829599522694E-01,
                      -2.068479571871339E-01,
                      -4.777798825282107E-02) | units.AU
  mercury.velocity = (9.121685164224357E-03,
                      2.532046488259770E-02,
                      1.231999349286273E-03) | (units.AU/units.day)
  mercury.radius = 2440.0 | units.km
  
  venus = planets[1]
  venus.mass = 48.685e23 | units.kg
  venus.position = (5.829554745446001E-01,
                    4.290508478691337E-01,
                    -2.773923286024579E-02) | units.AU
  venus.velocity = (-1.210377003565980E-02,
                    1.616561583873175E-02,
                    9.201724371252968E-04) | (units.AU/units.day)
  venus.radius = 6051.8 | units.km
  
  earth = planets[2]
  earth.mass = 5.97219e24 | units.kg
  earth.position = (2.914170806309638E-01,
                    -9.762380901507695E-01,
                    -8.838036341545362E-05) | units.AU
  earth.velocity = (1.621194321564609E-02,
                    4.839618876688764E-03,
                    -8.288236310894908E-07) | (units.AU/units.day)
  earth.radius = 6371.01 | units.km
  
  mars = planets[3]
  mars.mass = 6.4185e23 | units.kg
  mars.position = (-6.771052080404090E-01,
                   -1.358079217192958E+00,
                   -1.186843692826959E-02) | units.AU
  mars.velocity = (1.304721335990521E-02,
                   -5.060040804946996E-03,
                   -4.263506287278560E-04) | (units.AU/units.day)
  mars.radius = 3389.9 | units.km
  
  jupiter = planets[4]
  jupiter.mass = 1898.13e24 | units.kg
  jupiter.position = (-2.659291255441077E+00,
                      4.536415692760404E+00,
                      4.058673324381942E-02) | units.AU
  jupiter.velocity = (-6.601223938431606E-03,
                      -3.460147408421039E-03,
                      1.620752232688015E-04) | (units.AU/units.day)
  jupiter.radius = 71492.0 | units.km
  
  saturn = planets[5]
  saturn.mass = 5.68319e26 | units.kg
  saturn.position = (-6.149519345675703E+00,
                     -7.775694595902171E+00,
                     3.799392627119177E-01) | units.AU
  saturn.velocity = (4.071035079488449E-03,
                     -3.476237592282034E-03,
                     -1.012547571724344E-04) | (units.AU/units.day)
  saturn.radius = 60268.0 | units.km
  
  uranus = planets[6]
  uranus.mass = 86.8103e24 | units.kg
  uranus.position = (1.948219896924911E+01,
                     4.611924353502811E+00,
                     -2.352696714787791E-01) | units.AU
  uranus.velocity = (-9.347014358359851E-04,
                     3.643979378812646E-03,
                     2.560171100139048E-05) | (units.AU/units.day)
  uranus.radius = 25559.0 | units.km
  
  neptune = planets[7]
  neptune.mass = 102.41e24 | units.kg
  neptune.position = (2.731034114396837E+01,
                      -1.235250519790972E+01,
                      -3.750181405673675E-01) | units.AU
  neptune.velocity = (1.272538140363074E-03,
                      2.879023382710963E-03,
                      -8.821000085244238E-05) | (units.AU/units.day)
  neptune.radius = 24766.0 | units.km
  
  # coordinates transformation
  # with respect to solar system barycenter --> with respect to the Sun
  sun.position -= sun.position
  sun.velocity -= sun.velocity
  planets.position -= sun.position
  planets.velocity -= sun.velocity
  
  return sun, planets

def integrate_solar_system(sun, planets, time_end=5.0|units.yr, n_steps=500):
  """
  evolve the system using kepler_orbiters
  """
  
  converter = nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
  
  planets_around_sun = Kepler(converter,channel_type="sockets")
  
  # central particle
  planets_around_sun.central_particle.add_particles(sun[0:1])
  
  # to set the central particle at the center of the coordinate system
  #planets_around_sun.central_particle.position = (0.0, 0.0, 0.0) | units.AU
  #planets_around_sun.central_particle.velocity = (0.0, 0.0, 0.0) | units.kms
  
  # orbiters
  planets_around_sun.orbiters.add_particles(planets)
  planets_around_sun.commit_particles()
  
  # to change the integration method
  #planets_around_sun.parameters.method = 1
  #print planets_around_sun.get_method()
  
  channel_from_planetets_to_framework = planets_around_sun.orbiters.new_channel_to(planets)
  channel_from_sun_to_framework = planets_around_sun.central_particle.new_channel_to(sun)
  
  positions_sun = quantities.AdaptingVectorQuantity()
  positions_planets = quantities.AdaptingVectorQuantity()
  
  dt = time_end / float(n_steps)
  time = 0.0 | units.yr
  
  print " ** evolving solar system for", time_end.in_(units.yr), ", with time-step of", dt.in_(units.yr)
  print "    this might take a while"
  while time<=time_end:
    #print "\t", time.in_(units.yr)
    planets_around_sun.evolve_model(time)
    channel_from_planetets_to_framework.copy()
    channel_from_sun_to_framework.copy()
    positions_sun.append(sun.position)
    positions_planets.append(planets.position)
    time += dt
  print " **"
  
  planets_around_sun.stop()
  
  return positions_sun, positions_planets

def plot_solar_system(positions_sun, positions_planets):
  """
  simple plot to get the XY and XZ plane
  """
  
  fig = pyplot.figure(figsize=(10,7))
  
  # XY plane
  pyplot.subplot(121)
  pyplot.scatter(positions_sun[0,0,0].value_in(units.AU), 
                positions_sun[0,0,1].value_in(units.AU))
  for i,planet_i in enumerate(positions_planets[0,:,0]):
    pyplot.plot(positions_planets[:,i,0].value_in(units.AU), 
              positions_planets[:,i,1].value_in(units.AU))
  pyplot.axis([-45,45,-45,45])
  pyplot.gca().set_aspect('equal')
  pyplot.xlabel('x [AU]')
  pyplot.ylabel('y [AU]')
  
  
  # XZ plane
  pyplot.subplot(122)
  pyplot.scatter(positions_sun[0,0,0].value_in(units.AU), 
                positions_sun[0,0,2].value_in(units.AU))
  for i,planet_i in enumerate(positions_planets[0,:,0]):
    pyplot.plot(positions_planets[:,i,0].value_in(units.AU), 
              positions_planets[:,i,2].value_in(units.AU))
  pyplot.axis([-45,45,-10,10])
  pyplot.gca().set_aspect('equal')
  pyplot.xlabel('x [AU]')
  pyplot.ylabel('z [AU]')
  
  pyplot.tight_layout()
  pyplot.show()
  
  return

if __name__ in ('__main__','__plot__'):
  
  # initial conditions
  sun, planets = sun_and_planets()
  
  # integrate the solar system
  positions_sun, positions_planets = integrate_solar_system(sun, planets, time_end=10.0|units.yr, n_steps=1000)
  
  # simple plot of the results
  plot_solar_system(positions_sun, positions_planets)
  
  
