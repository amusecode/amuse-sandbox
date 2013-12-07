"""
   Nbody integration of N particles with a Salpeter initial mass
   function between Mmin and Mmax and with stellar evolution with
   metalicity z.
"""
from amuse.lab import *
from optparse import OptionParser
from matplotlib import pyplot

def plot_starplanet(f1, plot_xz):

    sizes = 100
    colors = 'b'
    alpha = 0.05
    bodies = read_set_from_file(f1, 'hdf5')
    pyplot.ion()
    pyplot.figure(figsize=(12,12))
    for bi in bodies.history:
        print "N=", len(bi)#, bi[0].age
        xb = bi[0:2].x.value_in(units.AU)
        yb = bi[0:2].y.value_in(units.AU)
        zb = bi[0:2].z.value_in(units.AU)
        xd = bi[2:].x.value_in(units.AU)
        yd = bi[2:].y.value_in(units.AU)
        zd = bi[2:].z.value_in(units.AU)
        pyplot.xlim(-150, 650)
        pyplot.ylim(-400, 400)
        if plot_xz:
            pyplot.scatter(xb, zb, c='r', lw=2)
            pyplot.scatter(xd, zd, sizes, colors, edgecolors = "none", alpha = alpha)
        else:
            pyplot.scatter(xb,yb, c='r', lw=2)
            pyplot.scatter(xd, yd, sizes, colors, edgecolors = "none", alpha = alpha)
#        pyplot.text(0.0, 0.95, "t="+str(int(t))+"Myr")
        pyplot.xlabel('x [AU]')
        pyplot.ylabel('y [AU]')
        pyplot.draw()
        pyplot.cla()
    pyplot.show()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="f1",  default = "HD106906_data.hdf5",
                      help="first file [%default]")
    result.add_option("--xz", 
                      dest="plot_xz", action='store_true', default = False,
                      help="plot x-z plane [%default]")

    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.AU, units.yr], 
                          precision = 10, prefix = "", 
                          separator = " [", suffix = "]")

    o, arguments  = new_option_parser().parse_args()
    plot_starplanet(**o.__dict__)
