from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from .interface import ClustorvmsInterface
from .interface import Clustorvms

class ClustorvmsInterfaceTests(TestWithMPI):
    
    def xtest1(self):
        instance = ClustorvmsInterface()
        result,error = instance.echo_int(12)
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()
    
