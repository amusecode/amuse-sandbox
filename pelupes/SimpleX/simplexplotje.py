import pylab
import numpy
import os

r=numpy.sqrt((results['x']-0.5)**2+(results['y']-0.5)**2+(results['z']-0.5)**2)
xion=numpy.log10(1-results['xion'])

pylab.scatter(r,xion)
pylab.show()
