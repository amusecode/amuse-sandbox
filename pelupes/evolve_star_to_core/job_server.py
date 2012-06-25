from amuse.rfi.core import *
import cPickle as pickle
from amuse.rfi.channel import AsyncRequestsPool
import inspect

def dump_and_encode(x,):
  return pickle.dumps(x)
def decode_and_load(x,):
  return pickle.loads(x)

class CodeImplementation(object):
   def exec_(self,arg):
     try:
       exec(arg)
       return 0
     except:  
       print ex
       return -1
   def _func(self,f,argin,kwargin,argout):
     try:
       func=decode_and_load(f)
       arg=decode_and_load(argin)
       kwarg=decode_and_load(kwargin)
       result=func(*arg,**kwarg)
       argout.value=dump_and_encode(result)
       return 0
     except Exception as ex:
       print ex
       argout.value=dump_and_encode(ex)
       return -1

class CodeInterface(PythonCodeInterface):    
    def __init__(self, **options):
        PythonCodeInterface.__init__(self, CodeImplementation, **options)
    
    @legacy_function
    def _func():
        function = LegacyFunctionSpecification()
        function.addParameter('func', dtype='string', direction=function.IN)
        function.addParameter('argin', dtype='string', direction=function.IN)
        function.addParameter('kwargin', dtype='string', direction=function.IN)
        function.addParameter('argout', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def exec_():
        function = LegacyFunctionSpecification()
        function.addParameter('arg', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function

    def func(self,f,*args,**kwargs):
        result,err=self._func( dump_and_encode(f),
                               dump_and_encode(args),
                               dump_and_encode(kwargs) )
        return decode_and_load(result[0]),err

    def async_func(self,f,*args,**kwargs):
        request=self._func.async(dump_and_encode(f),
                                 dump_and_encode(args),
                                 dump_and_encode(kwargs) )
        def f(x):
          result,err=x()
          return decode_and_load(result[0]),err
        request.add_result_handler( f )
        return request

class Job(object):
    def __init__(self, f, args, kwargs):
      self.f=f
      self.args=args
      self.kwargs=kwargs
      self._result=None
      self.request=None
      self.err=None

class JobServer(object):
    def __init__(self,hosts,channel_type="mpi",preamble=None):
      self.hosts=hosts
      self.job_list=[]
      self.idle_codes=[]
      self.last_finished_job=None
      self.channel_type=channel_type
      print "connecting hosts",
      i=0
      for host in hosts:
        i+=1; print i,
        self.idle_codes.append(CodeInterface(channel_type=self.channel_type,hostname=host))
      print
      if preamble is not None:
        for code in self.idle_codes:
          code.exec_(preamble)          
      self.pool=AsyncRequestsPool()
      print "AMUSE JobServer launched with", len(self.idle_codes),"threads"
    
    def submit_job(self,f,args=(),kwargs={}):
      job=Job(f,args,kwargs)
      self.job_list.append( job)
      if len(self.idle_codes)>0: 
          self._add_job(self.job_list.pop(), self.idle_codes.pop())        
      return job

    def wait(self):
      if len(self.pool)==0:
        return False
      else:
        self.pool.wait()
        return True   

    def waitall(self):
      while len(self.pool)>0:
        self.pool.wait()   

    def _finalize_job(self,request,job,code):
      job.result,job.err=request.result()
      if len(self.job_list)>0:
        self._add_job( self.job_list.pop(), code)
      else:
        self.idle_codes.append(code)  
      self.last_finished_job=job
    
    def _add_job(self,job,code):
      job.request=code.async_func(job.f,*job.args,**job.kwargs)
      self.pool.add_request(job.request,self._finalize_job, [job,code])
    
    def __del__(self):
      self.waitall()
      for code in self.idle_codes:
        code.stop()
