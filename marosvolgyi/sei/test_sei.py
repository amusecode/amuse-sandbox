from amuse.test.amusetest import TestWithMPI
import os
import sys
import numpy
import math

from interface import SeiInterface

class TestSeinterface(TestWithMPI):

    def test0(self):
        instance = SeiInterface()
        instance.initialization()
        instance.set_state(1,0,0,0,0,0)
        instance.evolve(0.1)
        print instance.get_state()
        instance.stop()
    
