import plummer2d
import time
#import threading
from Queue import Empty
from multiprocessing import Process, Queue, Event
from multiprocessing.managers import BaseManager

from amuse.units import nbody_system
from amuse.units import quantities

from amuse.community.hermite0.interface import Hermite


from optparse import OptionParser
class QueueManager(BaseManager):
    pass

QueueManager.register('get_data_queue')


def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-p", "--port", 
        default = 6504,
        dest="server_port",
        help="port on which server is listening",
        type="int"
    )
    return result
    
def run_hermite(queue):
    print "starting... 1"
    
    scale =  1 | nbody_system.length
    
    model = plummer2d.new_plummer_model_2D(20)
    code = Hermite()
    code.particles.add_particles(model)
    end_time = 10 | nbody_system.time
    print "starting..."
    try:
        for t in quantities.linspace(0 * end_time, end_time, 1000):
            code.evolve_model(t)
            print t
            queue.put(
                [code.particles.x / scale , code.particles.y / scale]
            )
    finally:
        print "stopping..."
        code.stop()
        print "done"
       
if  __name__ == '__main__':
    
    options, arguments = new_option_parser().parse_args()
    
    manager= QueueManager(address=('localhost', options.server_port), authkey='plotter')
    manager.connect()
    queue = manager.get_data_queue()
    
    run_hermite(queue)
