from amuse.units import units
import numpy
from amuse.community import *


uc = ConvertBetweenGenericAndSiUnits(
                1. | units.kpc,   
                1.e10 | units.MSun,
                1. | units.kms)
                
                
a=numpy.array([1.e20,2.e20,3.e20],'float32') | units.m

print a**2                

b=numpy.array([1.,2.,3.],'float32') | units.kpc

print uc.to_si(b)**2
