
import cPickle

from matplotlib import pyplot

from amuse.units import units,constants
from amuse.units.quantities import VectorQuantity,zero


f=open("snapshots/data.pkl","rb")

data=cPickle.load(f)

f.close()

t=data['time'][1:]
l=data['luminosity'][1:]

t=VectorQuantity.new_from_scalar_quantities(*t)
l=VectorQuantity.new_from_scalar_quantities(*l)

pyplot.plot(t.value_in(units.Myr),l.value_in(units.LSun),lw=2.)
pyplot.show()
