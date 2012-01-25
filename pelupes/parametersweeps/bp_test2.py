from python_interface import CodeInterface
from amuse.units import units
from amuse.rfi.channel import AsyncRequestsPool
import numpy

from matplotlib import pyplot

ncpu=4

results=[]

eccentricities=list(numpy.arange(0,0.5,0.1))

N=3
M=3

eccentricities= (0.5*i/(N)  for i in range(N+1))
eccnaxes= ( ((0.5*i)/N, 0.5+(1.5*j)/N ) for i in range(N+1) for j in range(M+1))


pool=AsyncRequestsPool()

def finalize_job(request,ecc,code):
  print "done with", ecc
  result,err=request.result()
  results.append(result)
  code.stop()
  del code
  add_job()

def add_job():
  try:
    ecc_a=eccnaxes.next()
  except:
    ecc_a=None
  if ecc_a is not None:  
    print "adding:",ecc_a
    ecc,a=ecc_a
    code=CodeInterface()
    request=code.async_func(
      m1=0.6897 | units.MSun,m2=0.20255 | units.MSun,m_planet=0.333 | units.MJupiter,
      r1=0.6489 | units.RSun,r2=0.22623 | units.RSun,r_planet=0.754 | units.RJupiter,
      ecc_binary=0.15944,P_binary=41.08| units.day,ecc_planet=ecc,a_planet=a | units.AU,
      pangle_planet=0., tend=200.| units.yr) 
    pool.add_request( request, finalize_job, [ecc_a,code])

for i in range(ncpu):
  add_job()

while len(pool)>0:
  pool.wait()

print results
