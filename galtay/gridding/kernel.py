import sys
import scipy.integrate
import numpy as np



SEVEN_OVER_EIGHTY = 7.0 / 80
FORTY_OVER_SEVEN_PI = 40.0 / ( 7 * np.pi )
GAUS_CONST_2D = 1 - np.exp( -40.0/7 )
SQRT_PI_OVER_TWO = np.sqrt(np.pi) / 2


class Kernel2D:
    """ This is the spline kernel used in gadget.  It is modified from the W_4
    spline kernel in http://adsabs.harvard.edu/abs/1985A%26A...149..135M  
    eq. 21 so that the kernel goes to zero at r=h instead of r=2h and 
    normalized so that a [volume,surface,line] integral over the whole kernel 
    in [3,2,1] dimensions gives 1.
    """

    def W( self, r, h ):
        """ The value of the kernel with smoothing lenght h at radius r."""
        x  = r / h
        x2 = x * x
        x3 = x2 * x

        norm = FORTY_OVER_SEVEN_PI / (h*h)

        if x <= 0.5:
            W = 1 - 6 * x2 + 6 * x3
        else:   
            W = 2 * (1-x)*(1-x)*(1-x)
            
        return W * norm 


    def return_value_at_samples( self, a, b, h, n ):
        """ Return n samples of the kernel between a and b. """
        if a < 0.0 or b < 0.0:
            print "a or b < 0.0"
            sys.exit(1)
        if a > h or b > h:
            print "a or b > h"
            sys.exit(1)

        r = np.linspace( a, b, n )
        W = np.zeros(n)
        for i in xrange(0,n):
            W[i] = self.W( r[i], h ) 

        return (r,W)



    def Wr( self, r, h, rmult=None, rpwr=None ):
        """ The value of the kernel with smoothing lenght h at radius r.
        Optionally multiply the result by some factor/power of r. """
        x  = r / h
        x2 = x * x
        x3 = x2 * x

        norm = FORTY_OVER_SEVEN_PI / (h*h)

        if x <= 0.5:
            W = 1 - 6 * x2 + 6 * x3
        else:   
            W = 2 * (1-x)*(1-x)*(1-x)

        if rmult==None:
            rmult = 1.0
        if rpwr==None:
            rpwr = 1.0


        return W * norm * rmult * r**rpwr




    def integrate_r_limits( self, a, b, h, rmult=2.0*np.pi, rpwr=1.0 ):
        """ Return the integral between r limits. """
        if a < 0.0 or b < 0.0:
            print "a or b < 0.0"
            sys.exit(1)
        if a > h or b > h:
            print "a or b > h"
            sys.exit(1)

        S = scipy.integrate.quad( self.Wr, a, b, args=(h,rmult,rpwr) )
        return S








class GaussianKernel2D:
    """ This is the gaussian approximation to the kernel above. """

    def W( self, r, h ):
        """ The value of the kernel with smoothing lenght h at radius r. """
        var =  h * h * SEVEN_OVER_EIGHTY
        norm = FORTY_OVER_SEVEN_PI / (h*h) / GAUS_CONST_2D

        r2 = r*r
        W = np.exp( - r2 / (2*var) )
        return W * norm


    def return_value_at_samples( self, a, b, h, n ):
        """ Return n samples of the kernel between a and b. """
        if a < 0.0 or b < 0.0:
            print "a or b < 0.0"
            sys.exit(1)
        if a > h or b > h:
            print "a or b > h"
            sys.exit(1)

        r = np.linspace( a, b, n )
        W = np.zeros(n)
        for i in xrange(0,n):
            W[i] = self.W( r[i], h ) 

        return (r,W)




    def Wr( self, r, h, rmult=None, rpwr=None ):
        """ The value of the kernel with smoothing lenght h at radius r.
        Optionally multiply the result by some factor/power of r. """
        var =  h * h * SEVEN_OVER_EIGHTY
        norm = FORTY_OVER_SEVEN_PI / (h*h) / GAUS_CONST_2D

        r2 = r*r
        W = np.exp( - r2 / (2*var) )

        if rmult==None:
            rmult = 1.0
        if rpwr==None:
            rpwr = 1.0
        
        return W * norm * rmult * r**rpwr




    def integrate_r_limits( self, a, b, h, rmult=2.0*np.pi, rpwr=1.0 ):
        """ Return the integral between r limits. """
        if a < 0.0 or b < 0.0:
            print "a or b < 0.0"
            sys.exit(1)
        if a > h or b > h:
            print "a or b > h"
            sys.exit(1)

        S = scipy.integrate.quad( self.Wr, a, b, args=(h,rmult,rpwr) )
        return S






    def return_erf_integrals( self, lims, h, DEBUG=False ):
        """ Returns the 1-D error function integrals. """
        var =  h * h * SEVEN_OVER_EIGHTY
        sqv = np.sqrt(2*var)
        AA = 1.0 / sqv

        tops = scipy.special.erf( AA * lims[1:]   ) 
        bots = scipy.special.erf( AA * lims[0:-1] )
        out = (SQRT_PI_OVER_TWO * sqv) * (tops - bots)

        return out

