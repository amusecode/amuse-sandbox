from amuse.units import units,constants
from amuse.datamodel import Particles 
import numpy, uvblue, bluered, spectra
from matplotlib import pyplot
from optparse import OptionParser

def main(temperature, radius, mass, resolution, uvmatch, brmatch, verbose):
    star = Particles(1)
    star.mass = mass | units.MSun
    star.temperature = temperature | units.K
    star.radius = radius | units.RSun
    metallicity = 0.02
    if uvmatch==None:
            uvmatch=True
    if brmatch==None:
            brmatch=True 
    eta=2
    lambda_resolution=1000
    resolution=resolution

    sp=spectra.Spectrum(star, metallicity,eta, uvmatch, brmatch, lambda_resolution, resolution, verbose)

    wavelenth = sp.lamb.number
    flux = sp.flux.number


    pyplot.semilogy(wavelenth, flux)
    pyplot.savefig('./spectrum.png')
    pyplot.show()

    print "Done!"

def new_option_parser():
    result = OptionParser()
    result.add_option("-T", dest="temperature", type="float", 
		      default = 5700.0, help="temperature of the star in K")
    result.add_option("-R", dest="radius", type="float", 
		      default = 1.0, help="radius of the star in RSun")
    result.add_option("-M", dest="mass", type="float", 
		      default = 1.0, help="mass of the star in MSun")
    result.add_option("--resolution", dest="resolution", type="float", 
		      default = None, help="Spectrum resolution")
    result.add_option("--no-uv", dest="uvmatch", action="store_false", 
	                 help="don't force UVBLUE library")
    result.add_option("--no-br", dest="brmatch", action="store_false", 
		             help="don't force BLUERED library")
    result.add_option("-v", dest="verbose", action="store_true", 
		             help="verbose mode")

    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

