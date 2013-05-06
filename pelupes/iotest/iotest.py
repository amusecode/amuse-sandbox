from amuse.datamodel import Particles

from amuse.io import write_set_to_file

from amuse.ic.plummer import new_plummer_model

p=new_plummer_model(500000)


print "writing"
write_set_to_file(p,"test.dat","amuse")
print "done"
