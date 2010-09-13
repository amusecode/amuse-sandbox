import numpy

from amuse.support.units import units
from amuse.support.units import constants

def mdot_leitherer1992_variant(star):
  m_msun=star.mass.value_in(units.MSun)
  l_lsun=star.luminosity.value_in(units.LSun)
  temp=star.temperature.value_in(units.K)
  z_zsun=1.
  return (units.MSun/units.yr).new_quantity(
    10**(-11.35+2.45*numpy.log10(l_lsun) - 1.1*numpy.log10(m_msun) +\
    -1.31*numpy.log10(temp)+0.80*numpy.log10(z_zsun)) )

def mdot_leitherer1992v_reimers(star):
  m_msun=star.mass.value_in(units.MSun)
  l_lsun=star.luminosity.value_in(units.LSun)
  r_rsun=star.radius.value_in(units.RSun)
  temp=star.temperature.value_in(units.K).clip(1000,8.e4)
  z_zsun=1.
  mdot=(10**(-12.+2.45*numpy.log10(l_lsun) - 1.1*numpy.log10(m_msun) +\
      -1.3*numpy.log10(temp)+0.80*numpy.log10(z_zsun)))
  i=numpy.where(temp<8000.)[0]
  mdot[i]=4.e-13*l_lsun[i]*r_rsun[i]/m_msun[i]
  return (units.MSun/units.yr).new_quantity(mdot)


def v_terminal_teff(star):
  t4=numpy.log10(star.temperature.value_in(units.K))-4.
  t4=t4.clip(0.,1.)
  return (30 | units.km/units.s) + ((4000 | units.km/units.s)*t4)

def e_supernova(star):
  i=numpy.where( (star.stellar_type>=13 | units.stellar_type) &
    (star.stellar_type<=15 | units.stellar_type) )[0] 
  n=len(star)
  e=numpy.array([0.]*n)
  e[i]=1.e51
  if(n == 1 ):
    return (units.erg).new_quantity(e[0])
  else:
    return (units.erg).new_quantity(e)


def lmech(star):
  lm=0.5*mdot_leitherer1992v_reimers(star)*v_terminal_teff(star)**2
  if(len(star) == 1 ):
    return lm[0]
  else:
    return lm  
