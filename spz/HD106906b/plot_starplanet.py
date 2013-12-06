"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
from amuse.lab import *
from optparse import OptionParser
from matplotlib import pyplot



def plot_starplanet(f1):

    bodies = read_set_from_file(f1, 'hdf5')
    for bi in bodies.history:
        print "N=", len(bi), bi[0].age
        x = bi.x.value_in(units.AU)
        y = bi.y.value_in(units.AU)
        pyplot.scatter(x,y, c='r', lw=2)
#        pyplot.text(0.0, 0.95, "t="+str(int(t))+"Myr")
        pyplot.xlabel('x [AU]')
        pyplot.ylabel('y [AU]')
    pyplot.show()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="f1",  default = "starlanet.hdf5",
                      help="first file [%default]")

    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.AU, units.yr], 
                          precision = 10, prefix = "", 
                          separator = " [", suffix = "]")

    o, arguments  = new_option_parser().parse_args()
    plot_starplanet(**o.__dict__)
