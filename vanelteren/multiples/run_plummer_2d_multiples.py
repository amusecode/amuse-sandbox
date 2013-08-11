import plummer2d
import time
import numpy
#import threading
from Queue import Empty
from multiprocessing import Process, Queue, Event
from multiprocessing.managers import BaseManager

from amuse.units import nbody_system
from amuse.units import quantities

from amuse.community.hermite0.interface import Hermite
from amuse.community.kepler.interface import Kepler
from amuse.community.smalln.interface import SmallN
from amuse.couple import encounters

import logging

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
    result.add_option(
        "-s", "--seed", 
        default = -1,
        dest="seed",
        help="random seed to use",
        type="int"
    )
    return result
    
def run_hermite(queue, seed):
    print "starting... 1"
    
    if seed > 0:
        numpy.random.seed(seed)
    
    scale =  1 | nbody_system.length
    
    model = plummer2d.new_plummer_model_2D(20)
    
    kepler = Kepler()
    kepler.initialize_code()
    
    smalln = SmallN()
    smalln.parameters.timestep_parameter = 0.1
    smalln.parameters.cm_index = 2001
    
        
    encounter_code = encounters.HandleEncounter(
        kepler_code =  kepler,
        resolve_collision_code = smalln,
        interaction_over_code = None
    )
    
    gravity_code = Hermite()
    code = encounters.Multiples(
        gravity_code = gravity_code,
        handle_encounter_code = encounter_code
    )
    model.radius = 0.01 | nbody_system.length
    code.particles.add_particles(model)
    end_time = 10 | nbody_system.time
    code.commit_particles()
    def plot_func(particles):
        queue.put(
                [particles.x / scale , particles.y / scale, 'r']
        )
    print "starting..."
    try:
        for t in quantities.linspace(0 * end_time, end_time, 1000):
            code.plot_func = plot_func
            code.evolve_model(t)
            print t
            particles = code.all_singles
            queue.put(
                [particles.x / scale , particles.y / scale, 'b']
            )
    finally:
        print "stopping..."
        gravity_code.stop()
        kepler.stop()
        smalln.stop()
        print "done"
       
if  __name__ == '__main__':
    
    options, arguments = new_option_parser().parse_args()
    
    logging.basicConfig(level=logging.ERROR)
    encounters.LOG_ENERGY.setLevel(logging.DEBUG)
    manager= QueueManager(address=('localhost', options.server_port), authkey='plotter')
    manager.connect()
    queue = manager.get_data_queue()
    
    run_hermite(queue, options.seed)
