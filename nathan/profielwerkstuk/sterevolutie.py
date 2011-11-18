"""
Sterevolutie

Dit python script genereert diverse figuren om de evolutie van sterren (massa 
tussen 0.1 en 100 zonsmassa's) te bestuderen.
"""

import sys
import numpy
from optparse import OptionParser

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
from amuse.plot import plot, semilogx, semilogy, loglog, xlabel, ylabel

from amuse.units import units
from amuse.datamodel import Particles
from amuse.community.sse.interface import SSE

def evolueer_ster(ster_massa):
    ster = Particles(1)
    ster.mass = ster_massa
    
    ster_evolutie_code = SSE()
    ster_evolutie_code.particles.add_particles(ster)
    ster_in_code = ster_evolutie_code.particles[0]
    
    tijd = [0.0] | units.yr
    lichtkracht = ster_in_code.luminosity.as_vector_with_length(1)
    straal = ster_in_code.radius.as_vector_with_length(1)
    temperatuur = ster_in_code.temperature.as_vector_with_length(1)
    massa = ster_in_code.mass.as_vector_with_length(1)
    
    while ster_in_code.stellar_type < 10 | units.stellar_type:
        ster_evolutie_code.evolve_model()
        tijd.append(ster_evolutie_code.model_time)
        lichtkracht.append(ster_in_code.luminosity)
        straal.append(ster_in_code.radius)
        temperatuur.append(ster_in_code.temperature)
        massa.append(ster_in_code.mass)
        
    ster_evolutie_code.stop()
    
    return tijd, lichtkracht, straal, temperatuur, massa
    
def maak_figuur(tijd, lichtkracht, straal, temperatuur, massa, figuurnaam):
    figure = pyplot.figure(figsize = (8, 12))
    figure.suptitle("Levensloop van een ster van "
        "{0:.3f} zonsmassa".format(massa[0].value_in(units.MSun)), fontsize=12)
    subplot1 = figure.add_subplot(4, 1, 1)
    semilogy(tijd[1:], lichtkracht[1:])
    xlabel('tijd', visible=False)
    ylabel('lichtkracht')
    pyplot.margins(0.05)
    
    subplot2 = figure.add_subplot(4, 1, 2, sharex=subplot1)
    semilogy(tijd[1:], straal[1:])
    xlabel('tijd', visible=False)
    ylabel('straal')
    pyplot.margins(0.05)
    
    subplot3 = figure.add_subplot(4, 1, 3, sharex=subplot1)
    semilogy(tijd[1:], temperatuur[1:])
    xlabel('tijd', visible=False)
    ylabel('temperatuur')
    pyplot.margins(0.05)
    
    figure.add_subplot(4, 1, 4, sharex=subplot1)
    plot(tijd[1:], massa[1:])
    xlabel('tijd')
    ylabel('massa')
    pyplot.margins(0.05)
    pyplot.subplots_adjust(hspace=0.001)
    for sub in [subplot1, subplot2, subplot3]:
        pyplot.setp(sub.get_xticklabels(), visible=False)
    pyplot.savefig(figuurnaam)
    
def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-m",
        "--massa", 
        dest="massa",
        type="float",
        default = 1.0,
        help="Massa van de ster in zonsmassa's"
    )
    return result


if __name__ == "__main__":
    options, arguments  = new_option_parser().parse_args()
    
    print "Simuleer de evolutie van een ster met massa:", options.massa | units.MSun
    tijd, lichtkracht, straal, temperatuur, massa = evolueer_ster(options.massa | units.MSun)
    
    figuurnaam = "massa_" + ("%0.3e"%(options.massa)).replace("+0","").replace("+","") + ".png"
    maak_figuur(tijd, lichtkracht, straal, temperatuur, massa, figuurnaam)
