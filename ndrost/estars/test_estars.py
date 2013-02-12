from amuse.community import *
from amuse.test.amusetest import TestWithMPI
from amuse.lab import *

from .interface import eStarsInterface
from .interface import eStars

class eStarsInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = eStars(channel_type='sockets')
        instance.initialize_code()
	converter = nbody.nbody_to_si(1|units.parsec, 1|units.MSun)
        plummer=new_plummer_model(100, converter)
        plummer.radius = 0.1|units.parsec
        plummer.red = 1.0
        plummer.green = 1.0
        plummer.blue = 1.0
	plummer.type = 0
	instance.particles.add_particles(plummer)

        instance.stop()
	
    
