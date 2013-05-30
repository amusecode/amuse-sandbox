# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Radiative Transfer in AMUSE: example with Stellar Evolution
# 
# In this exercise you will learn to design a simple numerical experiment using a radiative transfer code, coupled to stellar evolution (in highly idealized fashion).

# <markdowncell>

# ## imports
# 
# first we import a radiative tranfer code and a stellar evolution code:

# <codecell>

from amuse.community.sphray.interface import SPHRay
from amuse.community.sse.interface import SeBa

# <markdowncell>

# next, a simple initial condition for the gas:

# <codecell>

from amuse.ext.molecular_cloud import new_ism_cube

# <markdowncell>

# then the AMUSE units and constants module, and the Particle object:

# <codecell>

from amuse.units import units,constants
from amuse.datamodel import Particles

# <markdowncell>

# take a moment to consider the following:
# 
# 1. what are the pro's and con's of importing all the components seperately versus the:
# 
#         from amuse.lab import *
# 
#    statement used in the tutorial?

# <markdowncell>

# ## Getting started
# 
# The above imports make available a number of classes, which you are probably unfamiliar with - there are a number of ways to get information about these
# (without reading the manual!):
#    
#     print help(new_ism_cube)
# for the community codes, it is always a good idea to readup on the respective code papers:
# 
#     code=SPHRay()
#     code.print_literature_references()
# or print out its parameters:
# 
#     print code.parameters
# or:
# 
#     help(code.parameters)
# 
# 1. play around with these way to get information

# <codecell>

pass

# <markdowncell>

# ## coupling rad transfer and stellar evolution
# 
# consider the code below (which is slightly modified from the code you have seen before in the tutorial 6):

# <codecell>

tend=5. | units.Myr

gas = new_ism_cube(5000, 1 | units.kpc, 0.01 | (units.amu / units.cm**3))
gas.h_smooth = 0.1 | units.kpc
gas.xion = 0.00
sources = Particles(1,position = [0,0,0] | units.parsec, luminosity = 1 | 1.e50 / units.s, SpcType=0)

rad = SPHRay()
rad.parameters.box_size = 2.0 | units.kpc
rad.parameters.number_of_rays= 10000 | units.Myr**-1

rad.gas_particles.add_particles(gas)
rad.src_particles.add_particles(sources)
rad.evolve_model(tend)
scatter(rad.gas_particles.position.lengths().value_in(units.kpc), rad.gas_particles.xion)
rad.stop()

# <markdowncell>

# ## including stellar evolution
# 
# we are going to include sources based on a stellar evolution for a 30 MSun star.
# 
# 1. copy and change the above code to set the ionizing luminosity to a value calculated from the initial model of the stellar evolution code using the function you wrote in the previous exercise. For this, you need to generate the source with a mass instead of a luminosity. Note that a particle in a stellar evolution code also has an attribute luminosity, but this is the total blolometric luminosity.

# <codecell>

pass

# <markdowncell>

# ## time dependence
# 
# next, we want to add time dependence. This will be done by simply alternating stellar evolution and radiative transfer, updating the luminoisty of the source as the star evolves. For this we need a way to update the luminosity in SPHRay, this can be done using a channel:
#     
#     channel_to_rad=sources.new_channel_to(rad.src_particles)
#     channel_to_rad.copy_attributes(["x","y","z","luminosity","SpcType"])
# 
# similarly, we can copy the attributes needed for the calculation of the luminosity from the stellar evolution code se:
# 
#     channel_from_se=se.particles.new_channel_to(sources)
#     channel_from_se.copy_attributes(["radius","temperature"])
# 
# 1. change the code to evolve the model taking substeps dt, and make a plot every dtplot
# 2. change the code to co-evolve the stellar evolution
# 3. change the code to update the ionizing luminosity (using the above channel) every dt

# <codecell>

pass

# <markdowncell>

# if you have time you can consider the following:
#     
# 1. change the code so that instead of a single star, it contains radiation from a population of stars.
#     

# <codecell>

from amuse.ic.salpeter import new_salpeter_mass_distribution

