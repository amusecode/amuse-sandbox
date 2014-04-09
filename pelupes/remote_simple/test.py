from distributed_amuse import local_and_paddegat,local_only

from simple import simpleInterface,Simple

l=local_and_paddegat()

s=Simple(channel_type="distributed", distributedInstance=l.instance,
   redirection="none",label="paddegat",name_of_the_worker=remote_worker)

print s.timestwo(5.)

raw_input()

s.stop()
