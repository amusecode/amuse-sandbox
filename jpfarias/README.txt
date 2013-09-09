The package spectra.py contains a series of tools to have
access to the UVBLUE and BLUERED database.
For details about the database please visit the
UVBLUE and BLUERED homepage's:

http://www.inaoep.mx/~modelos/bluered/bluered.html
http://www.inaoep.mx/~modelos/uvblue/uvblue.html

This package depends of the uvblue.py and bluered.py 
package's written by Taro Sato (ubutsu@gmail.com), 
and the scipy library.

The database directory must be specified in the variables:
uvdir=<<path to UVBLUE database>>
brdir=<<path to BLUERED database>>

Files must be compressed using GunZip and have the following 
structure (that should be the original format if are downloaded 
from the official page).

UVBLUE example:
uvblue_path/t05000g50m10k2.flx.gz

BLUERED example:
bluered_path/spectrum.t05250g20m30k2.gz

In booth cases:

txxxxx  : Temperature [K]
gxx     : logg * 10
mxx     : negative mettalicity * 10 
pxx     : positive or zero mettaliciy *10
kx      :


#########################################################

User's Manual:

The package spectra contains the Spectrum class that will obtain 
the closest spectrum matching in the library. Where:

sp=spectra.Spectrum(star, metallicity=0.0, 
                    uvmatch=True,
                    brmatch=True, 
                    eta=2, 
                    lambda_resolution=None, 
                    resolution=None,
                    verbose=False)

star                : Is an AMUSE object, containing the quantity parameters: temperature, radius and mass.

metallicity         : star metallicity.

uvmatch             : Force UVBLUE database to match. It means that will choose the closest spectrum that
                      matches with the star.

brmatch             : Force BLUERED database to match. 
                      When uvmacth and brmatch are True, the spectrum will be the closest one that
                      match both libraries at the same time.
                      If uvmatch or brmatch are False (on at time), the spectrum will be the closest one
                      in the forced lybrary, and the rest will be approximated as a blackbody spectrum.  

lambda_resolution   : Lambda resolution where  R_l=lambda/d_lambda

resolution          : Resolution for the spectrum where : R_flux= lambda/d_lambda 
                      If resolution are not specified, will return the best resolution possible for
                      all the wavelength range, that is:

                      UVBlUE:  lambda_resolution=200000
                               resolution=50000

                      BLUERED: lambda_resolution=40000
                               resolution=10000

verbose             : If True will give screen information about the libraries used.
                      Default value False

The instance of Spectrum will have the following parameters:

flux                : Vector quantity with the flux of the closest spectrum matching the star, scaled by 
                      a factor f=BB(Teff)/BB(Tref) 
                      where BB(Teff) is a blackbody with the temperature of the star.
                            BB(Tref) is a blackbody with the temperature of the reference spectrum.

lamb                : Vector quantity containing the corresponding wavelength. Spaces between two points
                      follow the relation d_lambda=lambda/lambda_resolution

lambda_resolution   : Resolution used for lambda (could not be the specified by the user if the 
                      specified resolution is greater than the provided by the library.

resolution          : Resolution used for flux (could not be the specified by the user if the 
                      specified resolution is greater than the provided by the library.

Teff                : Temperature of the star.

logg                : computed log(g) of the star. 

metallicity         : metallicity of the star.

Tref                : Temperature in the library.

metallicity_ref     : metallicity in the library.

logg_ref            : log(g) in the library.


Some utility functions of the spectra package:

B_nu(nu,t)          : Returns the respective blackbody spectrum for a given frequency (nu)
                      at a given temperature (given by the Scalar Quantity t).

B_lambda(l,t)       : Returns the respective blackbody spectrum for a given wavelength (l)
                      at a given temperature (given by the Scalar Quantity t).

obtain_available_data() :  returns arrays containing the available spectrum in th UVBLUE and
                           BLUERED libraries. 
                           returns BLUERED_array, UVBLUE_array
                           whit format:
                           [temperature_in_Kelvin  logg*10 metallicity*10 k]

obtain_avaliable_metallicity() : Returns the available metallicity in the UVBLUE and BLUERED library.
                                 returns: BLUERED_met, UVBLUE_met

plot_available_data(metallicity=0.0):
                   
                   Plot the available data in a temperature-log(g) plane for a given metallicity.


