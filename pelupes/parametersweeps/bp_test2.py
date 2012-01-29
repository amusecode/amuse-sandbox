from python_interface import CodeInterface
from amuse.units import units
from amuse.rfi.channel import AsyncRequestsPool
import numpy

from matplotlib import pyplot

ncpu=20
allhosts=[ ("paddegat",4),
           ("koppoel",4),
           ("gaasp",4),
           ("biesbosch",4),
           ("bullewijk",4),
           ]

hosts=reduce(lambda x,y: x+[y[0]]*y[1],allhosts,[])
print hosts
#raise Exception

results=[]

N=10
M=10

amin=0.5
amax=3.

eccnaxes= ( ((0.5*i)/N, amin+((amax-amin)*j)/N ) for i in range(N+1) for j in range(M+1))


pool=AsyncRequestsPool()

def finalize_job(request,ecc,code):
  print "done with", ecc
  result,err=request.result()
  results.append(result)
#  code.stop()
#  del code
  add_job(code)

def add_job(code=None):
  try:
    ecc_a=eccnaxes.next()
  except:
    ecc_a=None
  if ecc_a is not None:  
    print "adding:",ecc_a
    ecc,a=ecc_a
    if code is None:
      code=CodeInterface(hostname=hosts.pop())
    request=code.async_func(
      m1=0.6897 | units.MSun,m2=0.20255 | units.MSun,m_planet=0.333 | units.MJupiter,
      r1=0.6489 | units.RSun,r2=0.22623 | units.RSun,r_planet=0.754 | units.RJupiter,
      ecc_binary=0.15944,P_binary=41.08| units.day,ecc_planet=ecc,a_planet=a | units.AU,
      pangle_planet=0., tend=100000.| units.yr) 
    pool.add_request( request, finalize_job, [ecc_a,code])

for i in range(ncpu):
  add_job()

while len(pool)>0:
  print "waiting for job to finish.."
  pool.wait()

f=open('results','w')
import cPickle
cPickle.dump(results,f)
f.close()
