"""
Reproduction of the Fig. 23 from Toomre & Toomre (1972, TT72)
  
  -- produces plot tt72_fig23.png
"""

import numpy

from amuse.units import units, constants, nbody_system
from amuse.datamodel import Particles, ParticlesSuperset
from amuse.community.kepler.interface import Kepler as Kepler_twobody
from amuse.units.optparse import OptionParser

from disk_flyby import integrate_2disks_flyby

from plot_disks import plot_all_snaps

pipi = 2.*numpy.pi
pi_180 = numpy.pi/180.

def get_tt72_disk(m=10.e11|units.MSun,
                  r_min=25.|units.kpc,
                  n_rings=[12,15,18,21,24,27,30,33,36,39,42,45],
                  r_rings_rel=[0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75],
                  disk_id='a',
                  eps=0.|units.m):
  """
  initialize disk ala TT72 (see the first paragraphs of sec. II and III of the paper)
  positions and velocities with respect to the central particle of mass m
  """
  disk = Particles()
  
  for i,ri in enumerate(r_rings_rel):
    
    disk_rad_i = Particles(n_rings[i])
    
    a = ri*r_min
    phi_i = numpy.linspace(0., pipi, num=n_rings[i], endpoint=False)
    
    disk_rad_i.x = a * numpy.cos(phi_i)
    disk_rad_i.y = a * numpy.sin(phi_i)
    disk_rad_i.z = 0. * a
    
    x_r = disk_rad_i.x/a
    y_r = disk_rad_i.y/a
    
    #vc = (constants.G*m/a)**0.5
    vc = ( constants.G*m*a**2/(a**2 + eps**2)**1.5 )**0.5
    disk_rad_i.vx = -vc * y_r
    disk_rad_i.vy =  vc * x_r
    disk_rad_i.vz = 0.0 * vc
    
    disk.add_particles(disk_rad_i)
  
  # test particles
  disk.mass = 0.|units.MSun
  
  # identification of the disk
  disk.id = disk_id
  
  return disk
  
def get_inclined_disk(disk, 
                      incl_deg=0., 
                      omega_deg=0.,
                      lon_deg=0.):
  """
  rotate particles by inclination, longitude of ascending node, and argument of pericenter
  """
  
  incl = incl_deg * pi_180
  omega = omega_deg * pi_180
  longitude = lon_deg * pi_180
  
  # get rotation matrix
  a1 = ([numpy.cos(longitude), -numpy.sin(longitude), 0.0], 
        [numpy.sin(longitude), numpy.cos(longitude), 0.0], 
        [0.0, 0.0, 1.0])
  a2 = ([1.0, 0.0, 0.0], [0.0, numpy.cos(incl), -numpy.sin(incl)], [0.0, numpy.sin(incl), numpy.cos(incl)])
  a3 = ([numpy.cos(omega), -numpy.sin(omega), 0.0], [numpy.sin(omega), numpy.cos(omega), 0.0], [0.0, 0.0, 1.0])
  rot = numpy.dot(numpy.dot(a1,a2),a3)
  
  # reshape positions and velocity vectors
  pos = disk.position.value_in(units.kpc)
  vel = disk.velocity.value_in(units.kpc/units.Myr)
  
  # rotate
  pos_rot = (numpy.dot(rot, pos.T)).T | units.kpc
  vel_rot = (numpy.dot(rot, vel.T)).T | (units.kpc/units.Myr)
  
  inclined_disk = Particles(len(disk))
  inclined_disk.mass = disk.mass
  inclined_disk.position = pos_rot
  inclined_disk.velocity = vel_rot
  inclined_disk.id = disk.id
  
  return inclined_disk
  
def get_galaxies_in_orbit(m_a=10.e11|units.MSun,
                          m_b=10.e11|units.MSun,
                          ecc=0.5,
                          r_min=25.|units.kpc,
                          t_start=None):
  """
  binary galaxy with orbit of given parameters
  -- if ecc=>1, start at t_start (p625, t_start=-10=-10*100Myr)
  -- if ecc<1, start at apocenter (p644)
  """
  
  converter=nbody_system.nbody_to_si(m_a+m_b,1|units.kpc)
  semi = r_min/(1.-ecc)
  
  # relative position and velocity vectors at the pericenter using kepler
  kepler = Kepler_twobody(converter)
  kepler.initialize_code()
  kepler.initialize_from_elements(mass=(m_a+m_b), semi=semi, ecc=ecc, periastron=r_min) # at periastron
  
  # evolve back till initial position
  if ( ecc<1. ):
    kepler.return_to_apastron()
  else:
    kepler.transform_to_time(t_start)
  
  # get time of the orbit
  t_orbit = kepler.get_time()
  
  rl = kepler.get_separation_vector()
  r = [rl[0].value_in(units.AU), rl[1].value_in(units.AU), rl[2].value_in(units.AU)] | units.AU
  vl = kepler.get_velocity_vector()
  v = [vl[0].value_in(units.kms), vl[1].value_in(units.kms), vl[2].value_in(units.kms)] | units.kms
  
  kepler.stop()
  
  # assign particle atributes
  galaxies = Particles(2)
  galaxies[0].mass = m_a
  galaxies[0].position = (0,0,0) | units.AU
  galaxies[0].velocity = (0,0,0) | units.kms
  galaxies[1].mass = m_b
  galaxies[1].position = r
  galaxies[1].velocity = v
  
  # identification
  galaxies[0].id = 'a0'
  galaxies[1].id = 'b0'
  
  galaxies.move_to_center()
  
  return galaxies, t_orbit

def print_option_parser(parser_options):
  """ print options into file """
  file_par = parser_options.fout.replace(".hdf5", ".par")
  f_par = open(file_par, 'w')
  for oi in parser_options.__dict__:
    #print oi, parser_options.__dict__[oi]
    f_par.write(oi+' '+str(parser_options.__dict__[oi])+'\n')
  f_par.close()

def new_option_parser():
  result = OptionParser()
  result.add_option("-n", 
                    dest="n_steps", type="int", default = 10,
                    help="number of steps [%default]")
  result.add_option("--fout", 
                    dest="fout", default=None,
                    help="output file [%default]")
  result.add_option("--i_a",
                    dest="incl_a", type="float", default = -60.0,
                    help="inclination of the disk a [%default]")
  result.add_option("--i_b",
                    dest="incl_b", type="float", default = 60.0,
                    help="inclination of the disk b [%default]")
  result.add_option("--o_a",
                    dest="omega_a", type="float", default = 30.0,
                    help="argument of periastron of the disk a [%default]")
  result.add_option("--o_b",
                    dest="omega_b", type="float", default = 30.0,
                    help="argument of periastron of the disk b [%default]")
  result.add_option("--eta_nbody",
                    dest="eta_nbody", type="float", default=0.0006,
                    help="Huayno eta parameter (~timestep) for nbody [%default]")
  result.add_option("--eps_factor",
                    dest="eps_factor", type="float", default=0.2,
                    help="scale for softenning -- eps_factor*R_periastron [%default]")
  result.add_option("--ecc",
                    dest="ecc", type="float", default=0.5,
                    help="eccentricity [%default]")
  result.add_option("--r_min",
                    dest="r_min", type="float", default=25., unit=units.kpc,
                    help="eccentricity [%default]")
  return result

  
if __name__ in ('__main__', '__plot__'):
  
  o, arguments  = new_option_parser().parse_args()
  
  ### write options to file
  if o.fout is not None:
    print_option_parser(o)
  
  ### get galactic disks
  disk_a = get_tt72_disk(disk_id='a', eps=o.eps_factor*o.r_min)
  disk_b = get_tt72_disk(disk_id='b', eps=o.eps_factor*o.r_min)
  
  ### get galaxies in orbit
  galaxies, t_orbit = get_galaxies_in_orbit(ecc=o.ecc)
  print " ** time to apocenter for galaxies-binary =", -t_orbit.in_(units.Myr)
  
  ### setup the disks with respect to the galaxies
  # because we rotate the disk with respect to the orbital plane,
  #   we use longitude of ascending node instead of argument of pericenter in TT72,
  #   where the orbit is rotated with respect to the disk (see Fig. 6 in TT72)
  disk_a_incl = get_inclined_disk(disk_a, 
                                  incl_deg=o.incl_a, 
                                  lon_deg=o.omega_a)
  disk_b_incl = get_inclined_disk(disk_b, 
                                  incl_deg=o.incl_b, 
                                  lon_deg=o.omega_b)
  disk_a_incl.position += galaxies[0].position
  disk_a_incl.velocity += galaxies[0].velocity
  disk_b_incl.position += galaxies[1].position
  disk_b_incl.velocity += galaxies[1].velocity
  
  ### particle superset
  disks = ParticlesSuperset([disk_a_incl, disk_b_incl])
  
  ### parameters for the integration
  # t_end = -2.5*t_orbit
  t_end = -2.24 * t_orbit   # approximately corresponds to the time showed in Fig. 23
  file_redir="none"
  huayno_eps2 = (o.eps_factor*o.r_min)**2
  huayno_inttype_parameter = 12
  
  # number of test particles per subset in gravity_in_subsets
  nsub = 20
  
  print " ** t_end", t_end.in_(units.Myr)
  print " ** HUAYNO parameters: eps2 =", huayno_eps2.in_(units.kpc**2), "eta =", o.eta_nbody, "inttype =", huayno_inttype_parameter 
  
  ### integrate flyby
  system_after_flyby = integrate_2disks_flyby(galaxies, disks, t_end, o.n_steps, 
                                              o.fout, file_redir, 
                                              o.eta_nbody,
                                              huayno_eps2,
                                              huayno_inttype_parameter,
                                              nsub,
                                              verbose=True)
  ### plotting
  plot_all_snaps(system_after_flyby, './', 'tt72_fig23.png')
  
  ###
  ### time for you to leave
  ###