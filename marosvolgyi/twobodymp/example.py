from amuse.legacy import *

from interface import twobodympInterface
from interface import twobodymp
import Gnuplot, Gnuplot.funcutils
import numpy

if __name__ == '__main__':
    instance = twobodympInterface()
    instance.initialization(1000)
    print "Precision set to {0}".format( instance.get_precision().precision)
    result = instance.new_particle(1.0, 0.0,
                                   1.0, 0.1, -0.1,
                                   -0.1, 0.1, -0.2)
    print "Added particle number {0}".format(result.index_of_the_particle)
    data=[]

    for i in numpy.arange(0.01, 10.0, 0.05):
        instance.evolve_system(i)
        res = instance.get_position(1)
        data.append([res.x_, res.y_, res.z_])

    g = Gnuplot.Gnuplot(debug=1)
    g.splot(data)
    s=raw_input()

    instance.stop()
    
