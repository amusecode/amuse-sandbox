import threading
from distributed_amuse import local_only,paddegat_only
from amuse.community.bhtree.interface import BHTree
from amuse.community.twobody.interface import TwoBody
from simple import simpleInterface,Simple

import time

def start():
  return simpleInterface(channel_type="distributed",copy_worker_code=True)
#  return BHTree(channel_type="distributed")

d=paddegat_only()
d.instance.use_for_distributed_workers()

print "starting"
raw_input()

s=start()

i=0
while True:
  #~ r=s.timestwo.async(i)
  #~ r.wait()
  #~ r=r.result()
  r=s.timestwo(i)
  #~ r=s.parameters.epsilon_squared
  print i
  i+=1
