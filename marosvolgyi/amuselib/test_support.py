from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from .interface import supportInterface

class supportInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = supportInterface()
        print instance.add([1,1,1,1,1],[1,1,1,1,1]).sum
        instance.stop()

    
    
