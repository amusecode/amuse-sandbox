import numpy

from matplotlib import pyplot

from amuse.io import read_set_from_file
from amuse.units import units,constants



if __name__=="__main__":

    pyplot.ion()
    pyplot.show()

    f=pyplot.figure(figsize=(8,6))

      
    directory="./snapshots"
    
    for i in range(91):
    
      gas=read_set_from_file(directory+"/gas-%6.6i"%i,"amuse")
      sink=read_set_from_file(directory+"/sink-%6.6i"%i,"amuse")
      if len(sink)>0:
        print len(gas),len(sink),sink.total_mass()+gas.total_mass()
      else:
        print len(gas),len(sink),gas.total_mass()
      gas.temperature=((2.2| units.amu)*gas.u/constants.kB).in_(units.K)
      pyplot.clf()
      pyplot.loglog(gas.density.value_in(units.amu/units.cm**3), gas.temperature.number,'r.')
      pyplot.xlim(1,1.e9)
      pyplot.ylim(1.,10000.)
      pyplot.xlabel("density (cm**-3)")
      pyplot.ylabel("T (K)")
      pyplot.draw()
#      raw_input()
