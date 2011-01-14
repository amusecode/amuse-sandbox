from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from .interface import twobodympInterface
from .interface import twobodymp

class twobodympInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = twobodympInterface()
        instance.initialization(100)
        instance.set_pos(1.0, 0.1, -0.1)
        instance.set_vel(-0.1, 0.1, -0.2)
        instance.setmu(1.0)
        print instance.get_pos()
        instance.evolve_system(0.1)
        print instance.get_pos()
        instance.stop()
    
