import plummer2d
import time
#import threading
from Queue import Empty
from multiprocessing import Process, Queue, Event

from amuse.units import nbody_system
from amuse.units import quantities

from amuse.community.hermite0.interface import Hermite

from matplotlib import pyplot
from matplotlib import animation

def run_hermite(queue, need_to_stop_event, finished_event):
    print "starting... 1"
    model = plummer2d.new_plummer_model_2D(20)
    code = Hermite()
    code.particles.add_particles(model)
    end_time = 10 | nbody_system.time
    print "starting..."
    for t in quantities.linspace(0 * end_time, end_time, 1000):
        code.evolve_model(t)
        print t
        queue.put(
            [code.particles.x , code.particles.y]
        )
        if need_to_stop_event.is_set():
            print "it is set!"
            break
    code.stop()
    time.sleep(2)
    finished_event.set()
    
def update_plot(x, queue, path):
    if queue.empty():
        return
    x = None
    y = None
    try:
        q = queue.get(False)
        if not q is None:
            x, y = q
    except Empty:
        pass
    if x is None:
        return
    path[0].remove()
    path[0] = plot.scatter(
        x / scale,
        y / scale
    )    
if  __name__ == '__main__':
    queue = Queue()
    stop_event = Event()
    finished_event = Event()
    
    background = Process(target=run_hermite, args = [queue, stop_event, finished_event])
    background.daemon = False
    background.start()
    
    scale =  1 | nbody_system.length
    figure = pyplot.figure()
    plot = figure.add_subplot(111)
    plot.set_xlim(-1,1)
    plot.set_ylim(-1,1)
    x0, y0 = queue.get()
    scatter = plot.scatter(
        x0 / scale,
        y0 / scale
    )    
    
    animation = animation.FuncAnimation(figure, update_plot, fargs = [queue,[scatter]], frames=100, interval = 50)
    try:
        pyplot.show()
    finally:
        stop_event.set()
        finished_event.wait()
        background.join()