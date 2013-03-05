from amuse.community import *
from amuse.test.amusetest import TestWithMPI
from amuse.lab import *

from .interface import eStarsInterface
from .interface import eStars

class eStarsTests(TestWithMPI):
    
    def test1(self):
        instance = eStars(channel_type='sockets')
        instance.initialize_code()
	converter = nbody.nbody_to_si(1|units.parsec, 1|units.MSun)
        plummer=new_plummer_model(10, converter)
        plummer.radius = 0.1|units.parsec
        plummer.red = 1.0
        plummer.green = 0.0
        plummer.blue = 0.0
        plummer.alpha = 1.0
	plummer.type = 1
	instance.particles.add_particles(plummer)

	instance.store_view(1|units.Myr)
        instance.stop()
	
    
