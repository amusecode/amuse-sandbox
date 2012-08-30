""" Plot test 1 from http://arxiv.org/abs/astro-ph/0603199.  
Note that the box is setup with the origin at the center.  """

import os
import glob
import numpy
import argparse
from matplotlib import pyplot

from amuse.support.units import units
from amuse.support.units import constants
from amuse.ext.radial_profile import radial_profile
from amuse.support.io import read_set_from_file

import gridding.gridder as G



# test1 specific variables
#========================================
class Iliev1Vars():

  def __init__( self ):
    self.Lbox = 13.2 | units.kpc
    self.RC = 2.59e-13 | (units.cm**3/units.s)
    self.nH = 1.0e-3 | (units.cm**(-3))
    self.Np = 5.0e48 | (units.s**(-1))
    self.trec = 1 / (self.RC * self.nH)
    self.xHII = 1.2e-3
    self.Rs = ( ( 3 * self.Np ) / \
                  ( 4 * numpy.pi * self.RC * self.nH**2 ) )**(1./3)

  def show( self ): 
    print
    print ' Iliev Test 1 parameters '
    print '--------------------------------------'
    print 'Lbox:      ', self.Lbox.in_(units.kpc)
    print 'Rstr:      ', self.Rs.in_(units.kpc)
    print 't_rec:     ', self.trec.in_(units.Myr)
    print 'nH:        ', self.nH
    print 'RC:        ', self.RC
    print 'xHII_init: ', self.xHII
    print


  def Rt( self, t ):
    """ Radius as a function of time. """
    t_ratio = t.value_in( units.Myr ) / \
        self.trec.value_in( units.Myr )
    Ri = self.Rs * ( 1 - numpy.exp( -t_ratio ) )**(1./3)
    return Ri


  def Vt( self, t ):
    """ Velocity as a function of time. """
    t_ratio = t.value_in( units.Myr ) / \
        self.trec.value_in( units.Myr )
    texp = numpy.exp( -t_ratio )  
    Vi = self.Rs / (3*self.trec) * texp / (1-texp)**(2./3)
    return Vi



CC = Iliev1Vars()



# parse command line
#------------------------------------------------------
def parse_command_line():

  parser = argparse.ArgumentParser(description='Plot test Iliev1.')
  
  parser.add_argument( 'snap_file',  
                       help='output file to plot' )
  
  parser.add_argument( 'plot_type',  
                       help='type of plot to make',
                       choices=['image','profile','quick_check'])
  
  args = parser.parse_args()

  return args


# safe file read
#------------------------------------------------------
def read_set_from_file_with_close(filename, format = 'csv', 
                                  **format_specific_keyword_arguments):
            import amuse.io.store
            processor = amuse.io.store.StoreHDF(filename, 
                                                False, 
                                                open_for_writing=False)
            r = processor.load()
            processor.close()
            return r



# plots averaged radial profiles
#------------------------------------------------------
def aplot(i, tag, xyc,xlim=None,ylim=None,ylabel=""):
  pyplot.figure(figsize=(7,7))
  for x,y,c in xyc:
    xa,ya=radial_profile(x,y,N=20)
    pyplot.semilogy(xa,ya,c)  
  pyplot.xlabel('L/Lbox')  
  pyplot.ylabel(ylabel)  
  if xlim is not None:
    pyplot.xlim(xlim)
  if ylim is not None:
    pyplot.ylim(ylim)
  pyplot.savefig(tag+'-%7.7i.png'%i)


# plots individual particle radial profiles
#------------------------------------------------------
def dplot(fnum, tag, xyc,xlim=None,ylim=None,ylabel=""):
  pyplot.figure(figsize=(7,7))
  for x,y,c in xyc:
    pyplot.semilogy( x, y, '.', color=c, markersize=1.0, ls='none' )
  pyplot.xlabel('L/Lbox')  
  pyplot.ylabel(ylabel)  
  if xlim is not None:
    pyplot.xlim(xlim)
  if ylim is not None:
    pyplot.ylim(ylim)
  pyplot.savefig('plots/' + tag + '-' + fnum + '.png' )



def frames( label='iliev1'):

  path = 'output/'

  lo = numpy.array( [-6.6, -6.6] )
  ng = 512
  dl = 13.6 / ng
  grid = G.Grid2D( lo, ng, dl )
  
  for infile in glob.glob( os.path.join(path, 'iliev1*') ):
    print "current file is: " + infile

    outfile = 'frames/' + infile.split('/')[1] + '.png'
    print 'outfile is: ' + outfile

    if os.path.isfile( outfile ): continue

    g = read_set_from_file_with_close( infile, 'amuse' )	

    G.kernel_gauss( grid, 
                    g.x.value_in(units.kpc), 
                    g.y.value_in(units.kpc),
                    g.h_smooth.value_in(units.kpc), 
                    g.mass.value_in(units.MSun) * g.xion, 
                    periodic=True, dtype=numpy.float64, clear=True )
 
    print 'min/max grid.dat: ', numpy.log10(grid.dat.min()), \
                                numpy.log10(grid.dat.max())

    f=pyplot.figure( figsize=(7.0,7.0) )
    pyplot.imshow( numpy.log10(grid.dat), vmin=-1.0, vmax=2.0 )
    pyplot.colorbar()
    pyplot.savefig( outfile )
    f.clear()
    pyplot.close(f)




def profile_plots( fname=None, g=None ):

  if fname==None and g==None:
    print 'Need either a file or some data to plot.'
    sys.exit(1)

  if fname and data:
    print 'Need one and only one.  Choose fname (plot from file),'
    print 'or g (plot passed in data)'
    sys.exit(1)

  if fname:
    print 'fname: ', fname
    fnum = fname.split( '-' )[-1]
    g = read_set_from_file(fname, 'amuse')
  else:
    # check g has the right attributes
  

  gamma = 5./3.
  mu = 1.| units.amu
  boxlen = 13.2 

  r = ( (g.x**2+g.y**2+g.z**2)**0.5 ).value_in(units.kpc)
  boxlen = 13.2
  r = r / (boxlen/2)
  
  print
  print 'x min/max: ', g.x.value_in(units.kpc).min(), \
                       g.x.value_in(units.kpc).max()
  print 'y min/max: ', g.y.value_in(units.kpc).min(), \
                       g.y.value_in(units.kpc).max()
  print 'z min/max: ', g.z.value_in(units.kpc).min(), \
                       g.z.value_in(units.kpc).max()
  print 'r min/max: ', r.min(), r.max()
  print
  print 'h min/max: ', g.h_smooth.value_in(units.kpc).min(), \
                       g.h_smooth.value_in(units.kpc).max() 
  print 'm min/max: ', g.mass.value_in(units.MSun).min(), \
                       g.mass.value_in(units.MSun).max() 


  T = ( (gamma-1)*g.u*mu/(1+g.xion)/constants.kB ).value_in(units.K)
  dens = (g.rho).value_in(units.amu/units.cm**3)
  pres = ( (gamma-1)*g.u*g.rho).value_in(units.g/units.cm/units.s**2 )


  print
  print 'xion min/max: ', g.xion.min(), g.xion.max()

  dplot(fnum, 'xion',
        ( (r, g.xion,   'r'),
          (r, 1-g.xion, 'g')  ),
        xlim=(0.,1.),
        ylim=(1.e-5,3.),
        ylabel="x, 1-x" )

  dplot(fnum, 'temp',
        ( (r,T,'r'),),
        xlim=(0.,1.),
        ylim=(1.0e3,1.e5),
        ylabel="T (K)" )

  dplot(fnum, 'pres',
        ( (r,pres,'r'),),
        xlim=(0.,1.),
        ylim=(1.e-17,1.e-14),
        ylabel="P (g/cm/s**2)" )

  dplot(fnum, 'rho',
        ((r,dens,'r'),),
        xlim=(0.,1.),
        ylim=(0.0001,0.01),
        ylabel="density (amu/cm**3)" )

  #dplot(i,'mach',((r/boxlen/2,mach,'r'),),
  #        xlim=(0.,1.),ylim=(1.e-5,10.),ylabel="mach number")

  #dplot(i,'vel',((r/boxlen/2,v,'r'),(r/boxlen/2,cs,'g'),),
  #        xlim=(0.,1.),ylim=(0.01,100.),ylabel="bulk, sound speed (km/s)")






def quick_check_data( t, mass, xion, verbose=True, debug=False ):

  """ Simply checks the ionized volume against the analytic solution. """

  r_ana = CC.Rt( t )

  # for this constant density field, ionized mass fraction is equal to 
  # the ionized volume fraction

  total_mass = numpy.sum( mass, dtype=numpy.float64 )
  HII_mass_init = total_mass * CC.xHII
  HII_mass_now = numpy.sum( mass * xion, dtype=numpy.float64 )

  ion_mass_frac = (HII_mass_now - HII_mass_init) / total_mass
  Vion = CC.Lbox**3 * ion_mass_frac
  r_Vion = ( ( 3 * Vion ) / ( 4 * numpy.pi ) )**(1./3)

  if verbose:
    print 
    print 't:          ', t.in_( units.Myr )
    print 'r Vion:     ', r_Vion.in_(units.kpc)
    print 'r analytic: ', r_ana.in_( units.kpc )
    print 'Vion/ana:   ' , (r_Vion/r_ana)
    print 
  
  if debug:
    print 'total mass: ', total_mass
    print 'HII mass init: ', HII_mass_init
    print 'HII mass now: ', HII_mass_now
    print 'ion mass frac: ', ion_mass_frac
    print 'ionized volume: ', Vion.in_(units.kpc**3)
    
  





def quick_check_file( fname, verbose=True ):

  """ Simply checks the ionized volume against the analytic solution. """

  fnum = fname.split( '-' )[-1]
  t = numpy.float( fnum ) / 1.0e3 | units.Myr  
  r_ana = Rt( t )

  g = read_set_from_file(fname, 'amuse')

  # for this constant density field, ionized mass fraction is equal to 
  # the ionized volume fraction

  total_mass = numpy.sum( g.mass, dtype=numpy.float64 )
  HII_mass_init = total_mass * CC.xHII
  HII_mass_now = numpy.sum( g.mass * g.xion, dtype=numpy.float64 )

  ion_mass_frac = (HII_mass_now - HII_mass_init) / total_mass
  Vion = CC.Lbox**3 * ion_mass_frac
  r_Vion = ( ( 3 * Vion ) / ( 4 * numpy.pi ) )**(1./3)

  if verbose:
    print ' Quick Check '
    print '--------------------------------------'
    print 'fname: ', fname
    print 't @ snap: ', t.in_( units.Myr )
    print 'total mass: ', total_mass
    print 'HII mass init: ', HII_mass_init
    print 'HII mass now: ', HII_mass_now
    print 'ion mass frac: ', ion_mass_frac
    print 'ionized volume: ', Vion.in_(units.kpc**3)
    print 'r Vion: ', r_Vion.in_(units.kpc)
    print 'r analytic: ', r_ana.in_( units.kpc )
    print 'Vion/ana:   ' , (r_Vion/r_ana)
    print 

  

  


if __name__=="__main__":

  args = parse_command_line()

  if args.plot_type == 'profile':
    profile_plots_file( args.snap_file )
  elif args.plot_type == 'quick_check':
    quick_check_file( args.snap_file )



