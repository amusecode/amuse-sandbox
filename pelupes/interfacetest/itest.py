from amuse.community.gadget2.interface import Gadget2

import logging
logging.basicConfig(level=logging.WARN)

g=Gadget2(channel_type="sockets")

print "flag:", g.parameters.no_gravity_flag

#g.stop()

g=Gadget2(channel_type="mpi")

print "flag2:", g.parameters.no_gravity_flag

g.stop()
