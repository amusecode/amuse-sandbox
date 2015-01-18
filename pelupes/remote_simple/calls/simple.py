import os

from amuse.community import *

from amuse.rfi.core import PythonCodeInterface

from amuse.support.interface import InCodeComponentImplementation

from amuse.units import units

import time

class simpleImplementation(object):
  
  def timestwo(self,xin,xout):
#   print "timestwo is running", os.getcwd()
#    time.sleep(10)
    xout.value=2*xin
    return 0

class simpleInterface(PythonCodeInterface):

    def __init__(self, **options):
        PythonCodeInterface.__init__(
            self,
            simpleImplementation,
            **options)  
  
    @legacy_function
    def timestwo():
        function = LegacyFunctionSpecification()
        function.addParameter('xin', dtype='float64', direction=function.IN,unit=units.none)
        function.addParameter('xout', dtype='float64', direction=function.OUT,unit=units.km)
        function.result_type = 'int32'
        return function
  
class Simple(InCodeComponentImplementation):
    def __init__(self, convert_nbody = None, **options):
        interface = simpleInterface(**options)
        
        InCodeComponentImplementation.__init__(
            self,
            interface,
            **options
        )
  
