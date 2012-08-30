""" 
Several type of gridding routines for learning and research applications. 
Particles are always projected along the z-axis.  Perform a rotation on the 
particle set to project along different directions. 
"""

import kernel
import numpy as np
import pylab as plt

from gridder_fc import fortran_code_module as fc


NG_MIN = 32   # minimum grid size to use
HSML_FAC = 3  # try to cover each hsml with this many cells



# Exceptions
#=======================================================================

class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class OffGridError(Error):
    """Exception raised for passing in particles that are off the grid. """
    def __init__(self, message ):
        Exception.__init__( self, message )








def nearest_grid_point(grid, xp, yp, w=None, clear=True, 
                       dtype=np.float32):
    """ Interpolate particle values to grid using NGP. """

    if not hasattr(grid, 'dat'):
        grid.dat = np.zeros( (grid.ng,grid.ng), dtype=dtype )

    if grid.dat.dtype != dtype:
        print 'warning dtype mismatch, clearing dat'
        grid.dat = np.zeros( (grid.ng,grid.ng), dtype=dtype )

    if clear:
        grid.dat[:] = 0.0

    if np.any( xp < grid.lo[0] ) or \
       np.any( xp > grid.hi[0] ) or \
       np.any( yp < grid.lo[1] ) or \
       np.any( yp > grid.hi[1] ):
           
        raise OffGridError, '\n Particle positions cant be off grid. \n '


    # just call histogram2d
    #-----------------------------------
    xbins = np.linspace( grid.lo[0], grid.hi[0], grid.ng+1 )
    ybins = np.linspace( grid.lo[1], grid.hi[1], grid.ng+1 )

    grid.dat += np.histogram2d( xp,
                                yp,
                                [xbins,ybins],
                                weights=w )[0]



def cloud_in_cell(grid, xp, yp, w=None, periodic=True, clear=True, 
                  dtype=np.float32):
    """ Interpolate particle values to grid using CIC. """
        
    if not hasattr(grid, 'dat'):
        grid.dat = np.zeros( (grid.ng,grid.ng), dtype=dtype )

    if grid.dat.dtype != dtype:
        print 'warning dtype mismatch, clearing dat'
        grid.dat = np.zeros( (grid.ng,grid.ng), dtype=dtype )

    if clear:
        grid.dat[:] = 0.0

    #if np.any( xp < grid.lo[0] ) or \
    #   np.any( xp > grid.hi[0] ) or \
    #   np.any( yp < grid.lo[1] ) or \
    #   np.any( yp > grid.hi[1] ):
    #       
    #    raise OffGridError, '\n Particle positions cant be off grid. \n '

    xbins = np.linspace( 0, grid.ng, grid.ng+1 )
    ybins = np.linspace( 0, grid.ng, grid.ng+1 )

    # positions in grid units.  
    #-----------------------------------
    xp_g = ( xp - grid.lo[0] ) / grid.dl
    yp_g = ( yp - grid.lo[1] ) / grid.dl

    # center of nearest grid cell 
    #-----------------------------------
    xp_ngp = np.trunc( xp_g ) + 0.5
    yp_ngp = np.trunc( yp_g ) + 0.5

    # distance between particle and ngp
    #-----------------------------------
    xp_dngp = ( xp_g - xp_ngp ) 
    yp_dngp = ( yp_g - yp_ngp ) 

    # add 4 2d histograms with different weights
    #-----------------------------------
    if w==None:
        w = 1.0

    # add to ngp
    #-------------------------------------
    Wij = ( 1.0 - np.abs(xp_dngp) ) * \
          ( 1.0 - np.abs(yp_dngp) )        
        
    grid.dat += np.histogram2d( xp_ngp,
                                yp_ngp,
                                [xbins,ybins], 
                                weights=Wij*w )[0]

    # ngp +/- x dir
    #-------------------------------------
    Wij = (       np.abs(xp_dngp) ) * \
          ( 1.0 - np.abs(yp_dngp) )        

    xpw = xp_ngp + np.sign(xp_dngp)
    if periodic:
        xpw[xpw>grid.ng-1] = xpw[xpw>grid.ng-1] - grid.ng
        xpw[xpw<0] = xpw[xpw<0] + grid.ng

    grid.dat += np.histogram2d( xpw, 
                                yp_ngp, 
                                [xbins,ybins], 
                                weights=Wij*w )[0]

    # ngp +/- y dir
    #-------------------------------------
    Wij = ( 1.0 - np.abs(xp_dngp) ) * \
          (       np.abs(yp_dngp) )        

    ypw = yp_ngp + np.sign(yp_dngp)
    if periodic:
        ypw[ypw>grid.ng-1] = ypw[ypw>grid.ng-1] - grid.ng
        ypw[ypw<0] = ypw[ypw<0] + grid.ng        
        
    grid.dat += np.histogram2d( xp_ngp, 
                                ypw,
                                [xbins,ybins], 
                                weights=Wij*w )[0]

    # ngp +/- x and +/- y dir
    #-------------------------------------
    Wij = (       np.abs(xp_dngp) ) * \
          (       np.abs(yp_dngp) )        

    grid.dat += np.histogram2d( xpw,
                                ypw,
                                [xbins,ybins], 
                                weights=Wij*w )[0]
    


def kernel_flat(grid, xp, yp, hsml, w=None, periodic=True, clear=True, 
                dtype=np.float32 ):
    """ Interpolate particle values to grid using a constant kernel. 
    Each particle is represented by a square with center at (xp[i],yp[i])
    and side length 2*hsml[i].  Each grid cell that contains any part of 
    the square gets an equal contribution from the particle.  
    """

    if not hasattr(grid, 'dat'):
        grid.dat = np.zeros( (grid.ng,grid.ng), dtype=dtype )

    if grid.dat.dtype != dtype:
        print 'warning dtype mismatch, clearing dat'
        grid.dat = np.zeros( (grid.ng,grid.ng), dtype=dtype )

    if clear:
        grid.dat[:] = 0.0

    if np.any( xp < grid.lo[0] ) or \
       np.any( xp > grid.hi[0] ) or \
       np.any( yp < grid.lo[1] ) or \
       np.any( yp > grid.hi[1] ):
           
        raise OffGridError, '\n Particle positions cant be off grid. \n '


    # number of particles to grid
    #-----------------------------------
    npar = hsml.size

    # positions in grid units.  
    #-----------------------------------
    xp_g = ( xp - grid.lo[0] ) / grid.dl
    yp_g = ( yp - grid.lo[1] ) / grid.dl

    # low/high corner of squares in grid units
    #-----------------------------------
    xsq_lo_g = ( xp_g - hsml / grid.dl ).astype(np.int32)
    xsq_hi_g = ( xp_g + hsml / grid.dl ).astype(np.int32)

    ysq_lo_g = ( yp_g - hsml / grid.dl ).astype(np.int32)
    ysq_hi_g = ( yp_g + hsml / grid.dl ).astype(np.int32)

    # number of cells covered in each dimension
    #-----------------------------------
    xsq_dl = ( xsq_hi_g - xsq_lo_g ) + 1 
    ysq_dl = ( ysq_hi_g - ysq_lo_g ) + 1 

    # total number of cells covered 
    #-----------------------------------
    prd_sq_dl = xsq_dl * ysq_dl

    #
    # [xy]sq_dl determines the extent of the particle
    # in each dimension in units of grid cells
    #   [xy]sq_dl     = # of cells covered
    #   [xy]sq_dl + 1 = # of cell edges
    #
    
    # put particles covering one cell on grid
    #-----------------------------------
    indx = np.where( prd_sq_dl == 1 )[0]
        
    if len(indx) != 0:
        xbins = np.linspace( 0, grid.ng, grid.ng+1 )
        ybins = np.linspace( 0, grid.ng, grid.ng+1 )
        
        grid.dat += np.histogram2d( xp_g[indx], 
                                    yp_g[indx], 
                                    [xbins,ybins], 
                                    weights=w )[0]


    # loop over particles (booooooooo)
    #-----------------------------------
    for ip in range(npar):

        if npar > 1000:
            if np.mod( ip,npar/10 ) == 0:
                print 'done ', ip, 'of ', npar

        # skip if in one cell
        #---------------------------------------------------
        if npar > 1:
            if prd_sq_dl[ip] == 1:
                continue

        # make a list of grid indices
        #---------------------------------------------------
        xi = np.array( [], dtype=np.int32 )
        yi = np.array( [], dtype=np.int32 )

        # loop over rows of region
        #---------------------------------------------------
        for iy in range( ysq_dl[ip] ):

            xi_tmp = xsq_lo_g[ip] + \
                np.arange( xsq_dl[ip], dtype=np.int32 )

            yi_tmp = ysq_lo_g[ip] + \
                iy * np.ones( xsq_dl[ip], dtype=np.int32 )

            xi = np.concatenate( (xi, xi_tmp) )
            yi = np.concatenate( (yi, yi_tmp) )
                
        # add particle to grid
        #---------------------------------------------------
        if periodic:
            xi[xi>grid.ng-1] = xi[xi>grid.ng-1] - grid.ng
            xi[xi<0] = xi[xi<0] + grid.ng

            yi[yi>grid.ng-1] = yi[yi>grid.ng-1] - grid.ng
            yi[yi<0] = yi[yi<0] + grid.ng

        # check number of indices = number of grid cells
        #---------------------------------------------------
        if xi.size != prd_sq_dl[ip]:
            print 'mismatch'
            sys.exit(1)

        # add to grid
        #---------------------------------------------------
        if w==None:
            grid.dat[xi,yi] += 1.0/prd_sq_dl[ip]
        else:
            grid.dat[xi,yi] += w[ip]/prd_sq_dl[ip]




            
def kernel_gauss(grid, xp, yp, hsml, w=None, periodic=True, clear=True, 
                 dtype=np.float32, use_fortran=True, verbose=False ):

    """ Interpolate particle values to grid using a gaussian kernel. 
    Each particle is represented by a square with center at (xp[i],yp[i])
    and side length 2*hsml[i].  Each grid cell that contains any part of 
    the square gets a gaussian weighted contribution.  
    """

    if not hasattr(grid, 'dat'):
        grid.dat = np.zeros( (grid.ng,grid.ng), dtype=dtype )

    if grid.dat.dtype != dtype:
        print 'warning dtype mismatch, clearing dat'
        grid.dat = np.zeros( (grid.ng,grid.ng), dtype=dtype )

    if clear:
        grid.dat[:] = 0.0


    # number of particles to grid
    #-----------------------------------
    npar = hsml.size

    # create kernel object
    #-----------------------------------
    kern = kernel.GaussianKernel2D()

    # from distribution of hsml and HSML_FAC
    # determine lowest/highest res grids we need
    #-----------------------------------
    hsml_max = hsml.max()
    dl_max = hsml_max / HSML_FAC
    ng_min_raw = grid.len / dl_max

    hsml_min = hsml.min()
    dl_min = hsml_min / HSML_FAC
    ng_max_raw = grid.len / dl_min

    # find first power of 2 greater than or equal to ng_min
    # and first power of 2 greater than or equal to ng_max
    #---------------------------------------------------------
    ng_min = 2
    while ng_min < ng_min_raw:
        ng_min = ng_min * 2

    ng_max = 2
    while ng_max < ng_max_raw:
        ng_max = ng_max * 2

    # dont let ng_min go below NG_MIN
    #---------------------------------------------------------
    if ng_min < NG_MIN: ng_min = NG_MIN

    if verbose:
        print
        print 'grid len:       ', grid.len
        print 'hsml min/max:   ', hsml_min, hsml_max
        print 'dl min/max:     ', dl_min, dl_max
        print 'ng_raw max/min: ', ng_max_raw, ng_min_raw
        print 'ng max/min:     ', ng_max, ng_min
        print 'ng_inp:         ', grid.ng
        print

    # create low res grids.  the grid
    # passed into the routine will be 
    # grid0
    #-----------------------------------
    grids = {'grid0':grid}
    ng = grid.ng
    dl = grid.dl
    igrid = 0
    while ng > ng_min:
        igrid = igrid + 1
        key = 'grid' + str(igrid)
        ng = ng / 2
        dl = dl * 2
        grids[key]     = Grid2D( grid.lo, ng, dl )
        grids[key].dat = np.zeros( (ng,ng), dtype=dtype )
        if verbose:
            print 'making low res grid: ', igrid, 'with ng = ', ng        

    ngrids = 1 + igrid
    if verbose:
        print 
        
    # positions in grid units.  
    #-----------------------------------
    xp_g = ( xp - grid.lo[0] ) / grid.dl
    yp_g = ( yp - grid.lo[1] ) / grid.dl

    # low/high corner of squares in grid units
    #-----------------------------------
    xsq_lo_g = ( xp_g - hsml / grid.dl ).astype(np.int32)
    xsq_hi_g = ( xp_g + hsml / grid.dl ).astype(np.int32)

    ysq_lo_g = ( yp_g - hsml / grid.dl ).astype(np.int32)
    ysq_hi_g = ( yp_g + hsml / grid.dl ).astype(np.int32)

    # number of cells covered in each dimension
    #-----------------------------------
    xsq_dl = ( xsq_hi_g - xsq_lo_g ) + 1 
    ysq_dl = ( ysq_hi_g - ysq_lo_g ) + 1 

    # total number of cells covered 
    #-----------------------------------
    prd_sq_dl = xsq_dl * ysq_dl

    #
    # [xy]sq_dl determines the extent of the particle
    # in each dimension in units of grid cells
    #   [xy]sq_dl     = # of cells covered
    #   [xy]sq_dl + 1 = # of cell edges
    #


    # to decide which grid to use.  we want use the hsml_min
    # and hsml_max attributes of the grids
    #--------------------------------------------------------
    ndone = 0
    for igrid in range(ngrids):

        key = 'grid' + str(igrid)
        
        # select particles for this grid
        #-----------------------------------------------
        if igrid == 0:
            c1 = hsml <  grids[key].hsml_max
            indx = np.where( c1 )
        elif igrid == ngrids-1:
            c1 = hsml >= grids[key].hsml_min
            indx = np.where( c1 )
        else:            
            c1 = hsml >= grids[key].hsml_min
            c2 = hsml <  grids[key].hsml_max
            indx = np.where( c1 & c2 )

        # report how many are being done
        #-----------------------------------------------
        npar_on_grid = indx[0].size
        if verbose:
            print '  doing ', npar_on_grid, 'particles on level', igrid

        # cycle if not particles on this grid
        #-----------------------------------------------
        if npar_on_grid == 0:
            continue

        # call fortran to do the heavy lifting
        #-----------------------------------------------
        if use_fortran:
            
            if w==None: w = np.ones(npar)
        
            if periodic:
                periodic_fc = 1
            else:
                periodic_fc = 0

            if verbose:
                verbose_fc = 1
            else:
                verbose_fc = 0
            
            # add particles to grid
            #-----------------------------------------------
            grids[key].dat[:] += fc.return_sph_grid( xp[indx], 
                                                     yp[indx], 
                                                     hsml[indx], 
                                                     w[indx],
                                                     grids[key].lo[0], 
                                                     grids[key].lo[1], 
                                                     grids[key].len, 
                                                     grids[key].ng,
                                                     periodic_fc,
                                                     verbose_fc )

        else:

            print '  you should really compile the f2py module!  '
            sys.exit(1)
            

        # add to running total
        #-----------------------------------------------
        ndone = ndone + len(indx[0])

    # report final total
    #-----------------------------------------------
    if verbose:
        print
        print 'ndone: ', ndone


    # now sum up low res grids onto high res grids
    # work in reverse from lowest to highest res
    #-----------------------------------------------
    for igrid in range(ngrids)[::-1]:
        ckey = 'grid' + str(igrid)
        fkey = 'grid' + str(igrid-1)

        if igrid == 0: 
            break
        else:
            grids[fkey].dat += fc.prolong_grid( grids[ckey].dat, 
                                                grids[fkey].ng )


        


    













class Grid2D:
    """ A two dimensional square grid patch. """

    def __init__( self, lo, ng, dl ):
        self.lo = lo
        self.ng = ng
        self.dl = dl
        self.len = ng * dl
        self.hi = self.lo + self.len

        self.hsml_min = dl * HSML_FAC
        self.hsml_max = 2 * dl * HSML_FAC


    def show_grid( grid, gfac=1.0, ax=None, show_log=True,
                   vmin=None, vmax=None, extent=None, cmap='jet',
                   interpolation='none', show_colorbar=True ):
        
        """ Makes an image of whatever is in self.dat. """

        # make axis if we dont have one
        #------------------------------------------------------
        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111, aspect='equal') 

        # set extent from grid
        #------------------------------------------------------
        if extent==None:
            extent=(grid.lo[0],grid.hi[0],grid.lo[1],grid.hi[1])

        # make image plot
        #------------------------------------------------------
        if show_log:

            im = ax.imshow( np.transpose( np.log10(grid.dat*gfac) ), 
                            origin='lower', 
                            extent=extent,
                            interpolation=interpolation, 
                            vmin=vmin,
                            vmax=vmax,
                            cmap=cmap )

        else:

            im = ax.imshow( np.transpose( grid.dat*gfac ), 
                            origin='lower', 
                            extent=extent,
                            interpolation=interpolation, 
                            vmin=vmin,
                            vmax=vmax,
                            cmap=cmap )

        # add colorbar
        #------------------------------------------------------        
        if show_colorbar:
            plt.colorbar(im)

        # return axis 
        #------------------------------------------------------        
        return ax




    def show_grid_lines(self, ax=None):

        # make axis if we dont have one
        #------------------------------------------------------
        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111, aspect='equal') 
            
        # add grid lines to the plot
        #------------------------------------------------------
        for ii in range(self.ng+1):
            ax.plot( [ self.lo[0] + ii*self.dl, self.lo[0] + ii*self.dl ],
                     [self.lo[1], self.hi[1]],
                     color='black' )
                       
            ax.plot( [ self.lo[0], self.hi[0] ], 
                     [ self.lo[1] + ii*self.dl, self.lo[1] + ii*self.dl ],
                     color='black' )

        # return axis 
        #------------------------------------------------------        
        return ax



    # Super function to plot halo grids 
    #=======================================================================
    def show_halo_grid( self, halo, gfac=1.0, c_grid=None, cfac=1.0,
                        show_fof_rad=True, show_sub_rads=True,
                        npar_in_sub_limit=100, pt_npar_limit=1,
                        pt_subrads=1, ax=None, show_log=True,
                        vmin=None, vmax=None, extent=None, cmap='jet',
                        interpolation=None, show_colorbar=True,
                        c_vmin=None, c_vmax=None, c_levels=None,
                        c_colors='black', fof_rad_color='black',
                        fof_rad_ls='dashed', sub_rad_color='red' ):
        
        """ Makes an image of whatever is in self.dat. """


        # make axis if we dont have one
        #------------------------------------------------------
        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111, aspect='equal') 

        # set extent from grid
        #------------------------------------------------------
        if extent==None:
            extent=( self.lo[0], self.hi[0],
                     self.lo[1], self.hi[1] )

        # make image plot
        #------------------------------------------------------
        if show_log:

            im = ax.imshow( np.transpose( np.log10(self.dat*gfac) ), 
                            origin='lower', 
                            extent=extent,
                            interpolation=interpolation, 
                            vmin=vmin,
                            vmax=vmax,
                            cmap=cmap )

        else:

            im = ax.imshow( np.transpose( self.dat*gfac ), 
                            origin='lower', 
                            extent=extent,
                            interpolation=interpolation, 
                            vmin=vmin,
                            vmax=vmax,
                            cmap=cmap )


        # add contours if we want
        #------------------------------------------------------        
        if c_grid:
        
            if show_log:

                cr = ax.contour( np.transpose( np.log10(c_grid.dat*cfac) ),
                                 origin='lower',
                                 extent=extent,
                                 vmin=c_vmin,
                                 vmax=c_vmax,
                                 levels=c_levels,
                                 colors=c_colors )
                
            else:

                cr = ax.contour( np.transpose( c_grid.dat*cfac ),
                                 origin='lower',
                                 extent=extent,
                                 vmin=c_vmin,
                                 vmax=c_vmax,
                                 levels=c_levels,
                                 colors=c_colors )
                
        # add radius for FOF halo
        #------------------------------------------------------
        if show_fof_rad:
            ax = plt.gca()
            cen = halo.fof.cen
            rad = halo.fof.rad
            circ = plt.Circle( (cen[0],cen[1]), radius=rad,
                               color=fof_rad_color, fill=False,
                               linestyle=fof_rad_ls )
            ax.add_patch(circ)

        # add radii for subhalos
        #------------------------------------------------------
        if show_sub_rads:
            ax = plt.gca()
            for i in range(halo.fof.nsubs):
                if halo.subs[i]['npar'][pt_npar_limit] >= npar_in_sub_limit:
                    cen = halo.subs[i]['cen']
                    rad = halo.subs[i]['rad'][pt_subrads]
                    circ = plt.Circle( (cen[0],cen[1]), radius=rad,
                                       color=sub_rad_color, fill=False )
                    ax.add_patch(circ)

        # add colorbar
        #------------------------------------------------------        
        if show_colorbar:
            plt.colorbar(im)


















