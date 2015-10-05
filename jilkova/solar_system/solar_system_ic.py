import numpy

from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particles, Particle
from amuse.community.kepler.interface import Kepler

pi_over_180 = numpy.pi/180.

def new_kepler():
  converter = nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
  kepler = Kepler(converter)
  kepler.initialize_code()
  return kepler

def get_position(mass_sun, mass_planet, ecc, semi, mean_anomaly, incl, argument, longitude, delta_t=0.|units.day):
  """
  cartesian position and velocity from orbital elements,
  where the orbit is evolved from given mean_anomaly 
  by time delta_t
  argument -- argument of perihelion
  longitude -- longitude of ascending node
  """
  kepler = new_kepler()
  kepler.initialize_from_elements(mass=(mass_sun+mass_planet),
                                  semi=semi,
                                  ecc=ecc,
                                  mean_anomaly=mean_anomaly)
  kepler.transform_to_time(time=delta_t)
  r = kepler.get_separation_vector()
  v = kepler.get_velocity_vector()
  
  kepler.stop()
  
  a1 = ([numpy.cos(longitude), -numpy.sin(longitude), 0.0], [numpy.sin(longitude), numpy.cos(longitude), 0.0], [0.0, 0.0, 1.0])
  a2 = ([1.0, 0.0, 0.0], [0.0, numpy.cos(incl), -numpy.sin(incl)], [0.0, numpy.sin(incl), numpy.cos(incl)])
  a3 = ([numpy.cos(argument), -numpy.sin(argument), 0.0], [numpy.sin(argument), numpy.cos(argument), 0.0], [0.0, 0.0, 1.0])
  A = numpy.dot(numpy.dot(a1,a2),a3)
  r_vec = numpy.dot(A,numpy.reshape(r,3,1))
  v_vec = numpy.dot(A,numpy.reshape(v,3,1))
  
  # for relative vectors
  r[0] = r_vec[0]
  r[1] = r_vec[1]
  r[2] = r_vec[2]
  v[0] = v_vec[0]
  v[1] = v_vec[1]
  v[2] = v_vec[2]
  
  return r,v

def get_sun_and_planets(delta_JD=0.|units.day):
  """
  eight planets of the Solar System
  as for JD = 2457099.500000000 = A.D. 2015-Mar-18 00:00:00.0000 (CT)
  http://ssd.jpl.nasa.gov/horizons.cgi
  """
  planets = Particles(8)
  
  # mass
  planets.mass = [3.302e23,
                  48.685e23,
                  5.97219e24,
                  6.4185e23,
                  1898.13e24,
                  5.68319e26,
                  86.8103e24,
                  102.41e24] | units.kg
  
  # eccentricity
  planets_ecc = [2.056263501026885E-01,
                 6.756759719005901E-03,
                 1.715483324953308E-02,
                 9.347121362500883E-02,
                 4.877287772914470E-02,
                 5.429934603664216E-02,
                 4.911406962716518E-02,
                 8.494660388602767E-03]
  
  # semi-major axis
  planets_semi = [3.870989725156447E-01,
                  7.233252880006816E-01,
                  1.000816989613834E+00,
                  1.523624142457679E+00,
                  5.203543088590996E+00,
                  9.547316304899041E+00,
                  1.915982879739036E+01,
                  2.997013749028780E+01] | units.AU
  
  # mean anomaly [degrees]
  planets_mean_anomaly = [2.256667460183225E+02,
                          3.096834722926926E+02,
                          6.970055236286768E+01,
                          5.013506750245609E+01,
                          1.213203242081277E+02,
                          1.423311616732398E+02,
                          2.079860620353052E+02,
                          2.712246916734600E+02]
  planets_mean_anomaly = numpy.array(planets_mean_anomaly) * pi_over_180
  
  # inclination [IN degrees]
  planets_inclination = [7.004026765179669E+00,
                         3.394480103844425E+00,
                         3.563477431351056E-03,
                         1.848403408106458E+00,
                         1.303457729562742E+00,
                         2.488017444885577E+00,
                         7.728000142736371E-01,
                         1.767720502209091E+00]
  planets_inclination = numpy.array(planets_inclination) * pi_over_180
  
  # Longitude of Ascending Node [OM degrees]
  planets_longitude = [4.831163083479358E+01,
                       7.663982595051040E+01,
                       1.775515437672556E+02,
                       4.951282677064384E+01,
                       1.005036717671826E+02,
                       1.135683875842263E+02,
                       7.388411509910506E+01,
                       1.317497218434830E+02]
  planets_longitude = numpy.array(planets_longitude) * pi_over_180
  
  # Argument of Perihelion [W degrees]
  planets_argument = [2.916964171964058E+01,
                      5.469102797401222E+01,
                      2.877495001117996E+02,
                      2.865420083537150E+02,
                      2.740725976811202E+02,
                      3.398666856578898E+02,
                      9.666856264946740E+01,
                      2.951871807292030E+02]
  planets_argument = numpy.array(planets_argument) * pi_over_180
  
  planets.name = ['Mercury',
                  'Venus',
                  'Earth',
                  'Mars',
                  'Jupiter',
                  'Satrun',
                  'Uranus',
                  'Neptune']
  
  ### to compare with JPL, mass of the Sun needs to be rescaled
  #mg_nasa = 1.32712440018e20 | (units.m**3 / units.s**2)
  #g_nasa = 6.67259e-11 | (units.m**3 / units.kg / units.s**2)
  #ms = mg_nasa / g_nasa
  
  sun = Particle()
  sun.name = 'Sun'
  #sun.mass = ms
  sun.mass = 1.0 | units.MSun
  sun.position = [0.,0.,0.] | units.AU
  sun.velocity = [0.,0.,0.] | units.kms
  
  # get the position and velocity vectors relative to sun 
  # by evolving in Kepler
  for i,ecc_i in enumerate(planets_ecc):
    r, v = get_position(sun.mass,
                        planets[i].mass,
                        planets_ecc[i],
                        planets_semi[i],
                        planets_mean_anomaly[i],
                        planets_inclination[i],
                        planets_longitude[i],
                        planets_argument[i],
                        delta_t=delta_JD)
    planets[i].position = r
    planets[i].velocity = v
    
  return sun, planets

def solar_system_in_time(time_JD=2457099.5|units.day):
  """
  Initial conditions of Solar system --
  particle set with the sun + eight planets,
  at the center-of-mass reference frame.

  Defined attributes: 
  name, mass, radius, x, y, z, vx, vy, vz
  """
  time_0 = 2457099.5 | units.day
  delta_JD = time_JD-time_0
  sun, planets = get_sun_and_planets(delta_JD=delta_JD)
  
  solar_system = Particles()
  solar_system.add_particle(sun)
  solar_system.add_particles(planets)
  
  ### to compare with JPL, relative positions and velocities need to be corrected for the
  # Sun's vectors with respect to the barycenter
  #r_s = (3.123390770608490E-03, -4.370830943817017E-04, -1.443425433116342E-04) | units.AU
  #v_s = (3.421633816761503E-06,  5.767414405893875E-06, -8.878039607570240E-08) | (units.AU / units.day)
  #print sun
  #print planets.position.in_(units.AU) + r_s
  #print planets.velocity.in_(units.AU/units.day) + v_s
  
  return solar_system

if __name__ in ('__main__', '__plot__'):
  
  solar_system = solar_system_in_time()
  print solar_system


