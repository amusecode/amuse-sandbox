import numpy

from amuse.units import units,constants

pi=numpy.pi
e=numpy.e
kB=constants.kB
h=constants.h
c=constants.c
Ry=constants.Rydberg_constant
sigma=constants.Stefan_hyphen_Boltzmann_constant

def B_nu(nu,t):
  return 2*h*nu**3/c**2 * 1./ (e**(h*nu/kB/t)-1)

# slow
def photon_flux(T,lowfreq=0.|units.s**-1,N=100000):
  nu=(numpy.arange(N+1)+1.)/N*(kB*T)/h*25.+ lowfreq
  n=pi*B_nu(nu,T)/(h*nu)
  return (n[1:]+n[:-1]).sum()/2*(nu[1]-nu[0])

def energy_flux(T,lowfreq=0.|units.s**-1,N=100000):
  nu=(numpy.arange(N+1)+1.)/N*(kB*T)/h*25.+ lowfreq
  b=pi*B_nu(nu,T)
  return (b[1:]+b[:-1]).sum()/2*(nu[1]-nu[0])


temp=map(lambda x: float( 27500 + x*2500 ), range(12)  ) # temperature [K]

tempmin=27500. 
tempmax=55000. 
dtemp=2500.
imintemp=0
imaxtemp=11

logg=map(lambda x: 3.+ x*0.25, range(8) )  # log g [cgs]

loggmin=3.
dlogg=0.25
iminlogg=0
imaxlogg=7

ionflux=numpy.array([[ 23.21      ,  22.93      ,  22.7       ,  22.52      ,
                       22.39      ,  22.28      ,  22.2       ,  22.13      ],
                     [ 23.66      ,  23.47      ,  23.3       ,  23.15      ,
                       23.01      ,  22.89      ,  22.8       ,  22.72      ],
                     [ 23.77613633,  23.82      ,  23.7       ,  23.6       ,
                       23.5       ,  23.41      ,  23.33      ,  23.26      ],
                     [ 23.97250147,  24.08      ,  23.98      ,  23.9       ,
                       23.84      ,  23.78      ,  23.72      ,  23.67      ],
                     [ 24.14654933,  24.14654933,  24.2       ,  24.13      ,
                       24.08      ,  24.04      ,  24.        ,  23.97      ],
                     [ 24.30232564,  24.30232564,  24.38      ,  24.32      ,
                       24.28      ,  24.24      ,  24.21      ,  24.19      ],
                     [ 24.44293602,  24.44293602,  24.44293602,  24.48      ,
                       24.44      ,  24.41      ,  24.39      ,  24.37      ],
                     [ 24.57080494,  24.57080494,  24.57080494,  24.62      ,
                       24.58      ,  24.56      ,  24.54      ,  24.52      ],
                     [ 24.68785355,  24.68785355,  24.68785355,  24.74      ,
                       24.71      ,  24.69      ,  24.67      ,  24.66      ],
                     [ 24.79562447,  24.79562447,  24.79562447,  24.79562447,
                       24.82      ,  24.8       ,  24.79      ,  24.78      ],
                     [ 24.89537119,  24.89537119,  24.89537119,  24.89537119,
                       24.92      ,  24.9       ,  24.89      ,  24.88      ],
                     [ 24.98812331,  24.98812331,  24.98812331,  24.98812331,
                       25.01      ,  24.99      ,  24.98      ,  24.97      ]])
# from  OSTAR2002 Model Atomspheres
# solar metallicity     
# Number of ionizing photons in H I Lyman continuum
#               (log number per sec per cm2 at stellar surface)
# missing values in orignal table are calculated from blackbody
# ionflux[temp, logg]

def interpolate_ionizing_flux( logg, t):
  if t < tempmin:
    return 0. | ( units.cm**-2 * units.s**-1 )
  if t >= tempmax:
    return photon_flux( t | units.K,lowfreq=(c*Ry))
  
  ilogg_=(logg-loggmin)/dlogg
  ilogg=numpy.long( numpy.floor(ilogg_) )
  dl=ilogg_-ilogg
  
  itemp_=(t-tempmin)/dtemp
  itemp=numpy.long( numpy.floor(itemp_) )
  dt=itemp_-itemp

  itemp1=itemp+1
  ilogg1=ilogg+1

  flux=ionflux[itemp1,ilogg1] * dl     * dt     + \
       ionflux[itemp1,ilogg]  * (1-dl) * dt     + \
       ionflux[itemp,ilogg1]  * dl     * (1-dt) + \
       ionflux[itemp,ilogg]   * (1-dl) * (1-dt)
  return 10**flux | ( units.cm**-2 * units.s**-1 )

def ionizing_flux( Mstar, Rstar, Teff):
  g=constants.G*Mstar/Rstar**2
  logg=numpy.log10( g.value_in(units.cm/units.s**2)  )
  t=Teff.value_in(units.K)
  return interpolate_ionizing_flux(logg,t)

def ionizing_luminosity(Mstar,Rstar,Teff):
  return 4*numpy.pi*Rstar**2*ionizing_flux(Mstar,Rstar,Teff)

if __name__=="__main__":
  Mstar=12.09 | units.MSun
  Rstar=4.6| units.RSun
  Teff=30000. | units.K
  print ionizing_flux(Mstar,Rstar,Teff)
  print numpy.log10( ionizing_flux(Mstar,Rstar,Teff).value_in(units.cm**-2 * units.s**-1 ))
  print ionizing_luminosity(Mstar,Rstar,Teff).value_in(units.s**-1)

  l=ionizing_luminosity(Mstar,Rstar,Teff)
  m=1000. | units.MSun
  
  print (m/(1|units.amu)/l).in_(units.Myr)
  
