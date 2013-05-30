# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # AMUSE units: the blackbody spectrum
# 
# this is a small exercise to familiarize yourself with AMUSE units usage, learn to write simple AMUSE utility functions.

# <markdowncell>

# ## (1) first some definitions:
# 
# 1. check that you understand the following definitions 
# 2. there is at least one bug below - if you dont find it don't worry about it just yet!
# 3. what are the differences or similarities between an AMUSE unit, quantity and constant?

# <codecell>

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
  return 2*h*nu**4/c**2 * 1./ (e**(h*nu/kB/t)-1)

def B_lambda(l,t):
  return 2*h*c**2/l**5*1./(e**(h*c/(l*kB*t))-1)

def energy(nu):
  return constants.h*nu
  
def freq(energy):
  return energy/constants.h

def freq_from_wavenumber(k):
  return c*k

def wavelength(nu):
  return c/nu

def freq_from_wavelength(l):
  return c/l
  
def wiens_lambda_max(T):
  b = 2897768.6 | units.nano(units.m)* units.K
  return b/T

def wiens_T_from_lambda_max(l):
  b = 2897768.6 | units.nano(units.m)* units.K
  return b/l
    
def total_bolometric_flux(T):
  return sigma*T**4  

# <markdowncell>

# ## (2) using the above, calculate the following:
# 
# 1. the wavelength of the maximum for a T=10000K blackbody spectrum
# 2. the bolometric flux for a T=10000K blackbody spectrum
# 3. the total luminosity of the sun (Teff = 5780 K), in units.LSun!

# <codecell>

pass

# <markdowncell>

# ## (3) energy flux calculation
# 
# 1. what does the energy_flux function below do?
# 2. what are the T, lowfreq, N arguments? 
# 3. calculate the flux for some values of T (like 5500 K, 10000K, 50000K), try to convert to W/m**2, check against total bolometric flux above,
# 4. (spoiler alert) if you had spotted the bug: repeat 3. using the original buggy version of B_nu (with nu^4 instead of nu^3)..
# 5. consider the numpy.trapz call: observe the use of .number and .unit, can you understand what this does?
# 6. make a copy of the function and try to remove the .number and .unit - does the function still work? Can you think of a reason to use the form given here? (hint: what are b and b.number?)

# <codecell>

def energy_flux(T,lowfreq=0.|units.s**-1,N=10000):
  nu=(numpy.arange(1,N+1))/N*(kB*T)/h*25.+ lowfreq
  b=pi*B_nu(nu,T)
  return numpy.trapz(b.number,x=nu.number)| (b.unit*nu.unit)

# <codecell>

pass

# <markdowncell>

# ## (4) the photon flux
# 
# 1. write a function that calculates the photon flux, with an optional lower frequency bound
# 2. what would the lower bound have to be to calculate the ionizing flux?

# <codecell>

def photon_flux():
    pass

# <markdowncell>

# ## (5) stellar evolution
# 
# 1. instantiate a stellar evolution code (e.g. SSE or SeBa), and generate a stellar model for a 30. MSun star
# 2. print out its attributes
# 3. write a function that calculates the luminosity in ionizing photons Lion for a star (assuming of course it radiates like a blackbody!)
# 4. make a plot of Lion vs time

# <codecell>

from amuse.community.sse.interface import SSE
from amuse.community.seba.interface import SeBa
from amuse.datamodel import Particle

