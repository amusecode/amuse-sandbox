from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from .interface import MiInterface

class MiInterfaceTests(TestWithMPI):
    
    def xxtest1(self):
        instance = MiInterface()
        instance.stop()
    
