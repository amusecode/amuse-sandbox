"""
  a simple remote worker using AMUSE communication channel

usage:

  remote=RemoteCodeInterface()

  runs a remote AMUSE python worker process which can execute arbitrary 
  code remotely.

  remote.assign(express, arg): assigns the local variable (any pickleable) 
                  to the remote expression express (string)
  remote.execute(express): execute the string express
  remote.evaluate(express): evaluate the string expression express and 
                  sends it back
  remote.func(self,f,*args,**kwargs): execute local function f remotely 
                  with args and kwargs (f, args and kwargs must be 
                  pickleable

"""

import cPickle as pickle
from amuse.rfi.core import *


def dump_and_encode(x):
  return pickle.dumps(x,0)
def decode_and_load(x):
  return pickle.loads(x.encode("latin-1"))

class RemoteCodeImplementation(object):
   def __init__(self):
     self.scope={}
     self.scope['dump_and_encode']=dump_and_encode
     self.scope['decode_and_load']=decode_and_load

   def _exec(self,express):
     try:
       exec express in self.scope
       return dump_and_encode(None)
     except Exception as ex:
       return dump_and_encode(ex)
   def _eval(self,express,argout):
     try:
       self.scope.update(dict(express=express))
       exec "argout="+express in self.scope
       argout.value=eval("dump_and_encode(argout)",self.scope)
       return dump_and_encode(None)
     except Exception as ex:
       argout.value=dump_and_encode("")
       return dump_and_encode(ex)
   def _assign(self,lhs,argin):
     try:
       self.scope.update(dict(argin=argin))
       exec lhs+"=decode_and_load(argin)" in self.scope
       return dump_and_encode(None)
     except Exception as ex:
       return dump_and_encode(ex)
   def _func(self,f,argin,kwargin,argout):
     try:
       self.scope.update(dict(f=f,argin=argin,kwargin=kwargin))
       exec "func=decode_and_load(f)" in self.scope
       exec "arg=decode_and_load(argin)" in self.scope
       exec "kwarg=decode_and_load(kwargin)" in self.scope
       exec "result=func(*arg,**kwarg)" in self.scope
       argout.value=eval("dump_and_encode(result)",self.scope)
       return dump_and_encode(None)
     except Exception as ex:
       argout.value=dump_and_encode("")
       return dump_and_encode(ex)

class RemoteCodeInterface(PythonCodeInterface):    
    def __init__(self, **options):
        PythonCodeInterface.__init__(self, RemoteCodeImplementation, **options)
    
    @legacy_function
    def _func():
        function = LegacyFunctionSpecification()
        function.addParameter('func', dtype='string', direction=function.IN)
        function.addParameter('argin', dtype='string', direction=function.IN)
        function.addParameter('kwargin', dtype='string', direction=function.IN)
        function.addParameter('argout', dtype='string', direction=function.OUT)
        function.result_type = 'string'
        return function

    @legacy_function
    def _exec():
        function = LegacyFunctionSpecification()
        function.addParameter('arg', dtype='string', direction=function.IN)
        function.result_type = 'string'
        return function

    @legacy_function
    def _eval():
        function = LegacyFunctionSpecification()
        function.addParameter('arg', dtype='string', direction=function.IN)
        function.addParameter('argout', dtype='string', direction=function.OUT)
        function.result_type = 'string'
        return function

    @legacy_function
    def _assign():
        function = LegacyFunctionSpecification()
        function.addParameter('lhs', dtype='string', direction=function.IN)
        function.addParameter('argin', dtype='string', direction=function.IN)
        function.result_type = 'string'
        return function

    def execute(self,express):
        err=decode_and_load( self._exec(express)[0] )
        if err:
          raise err

    def assign(self,lhs,arg):
        err=decode_and_load( self._assign(lhs, dump_and_encode(arg))[0] )
        if err:
          raise err

    def evaluate(self,express):
        result,err=self._eval(express)
        err=decode_and_load( err[0])
        if err :
          raise err
        return decode_and_load(result[0]) 

    def func(self,f,*args,**kwargs):
        result,err=self._func( dump_and_encode(f),
                               dump_and_encode(args),
                               dump_and_encode(kwargs) )
        err=decode_and_load( err[0])
        if err :
          raise err
        return decode_and_load(result[0])

    def async_func(self,f,*args,**kwargs):
        request=self._func.async(dump_and_encode(f),
                                 dump_and_encode(args),
                                 dump_and_encode(kwargs) )
        def f(x):
          result,err=x()
          err=decode_and_load( err[0])
          if err :
            raise err
          return decode_and_load(result[0])
        request.add_result_handler( f )
        return request
