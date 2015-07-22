from amuse.units import units

from matplotlib  import pyplot

import cPickle

import numpy

c="bgry"

filenames=["t_rms_0.00",
           "t_rms_0.20",
          "t_rms_1.00",
          "t_rms_5.00"]

labels=["0.xM_jup","0.2xM_jup","1.xM_jup","5.xM_jup"]

for i,filename in enumerate(filenames):
  f=open(filename,"r")
  t,rms=cPickle.load(f)
  f.close()
  
  pyplot.semilogy(t.value_in(units.yr),rms.value_in(units.AU),c[i],label=labels[i])

pyplot.xlim(0,10000)
pyplot.ylim(1.e-5,10.)

pyplot.xlabel("time (yr)")
pyplot.ylabel("rms position spread (AU)")

pyplot.legend(loc="upper left")

pyplot.savefig("t_rms.eps")  
pyplot.show()
