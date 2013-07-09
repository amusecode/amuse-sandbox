from Queue import Empty
from multiprocessing import Process, Queue, Event
from multiprocessing.managers import BaseManager

from matplotlib import pyplot
from matplotlib import animation

from optparse import OptionParser

class QueueManager(BaseManager):
    pass
    
def run_manager(queue,  stop_event, port, authkey = 'plotter'):
    manager = QueueManager(address=('', port), authkey=authkey)
    manager.register('get_data_queue', callable = lambda : queue)
    manager.register('shutdown', callable = lambda: server.shutdown())
    server = manager.get_server()
    server.serve_forever()
    
class QueuedPlotter(object):
    
    def __init__(self, queue, port, authkey = 'plotter'):
        self.data_queue = queue
    
    def get_data_queue(self):
        return self.data_queue
        
def update_plot(x, queue, path):
    if queue.empty():
        return
    x = None
    y = None
    c = None
    try:
        q = queue.get(False)
        if not q is None:
            x, y, c = q
    except Empty:
        pass
    if x is None:
        return
    
    if c == 'r':
        plot2.scatter(x, y, c = c)
    else:
        path[0].remove()
        path[0] = plot1.scatter(
            x,
            y,
            c = c
        )   
        

def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-p", "--port", 
        default = 6504,
        dest="server_port",
        help="port to listen on for events",
        type="int"
    )
    return result
    
if  __name__ == '__main__':
    
    options, arguments = new_option_parser().parse_args()
    
    stop_event = Event()
    
    queue = Queue()
    
    background = Process(target=run_manager, args = [queue, stop_event, options.server_port])
    background.daemon = True
    background.start()
    
    figure = pyplot.figure(figsize=(5,5))
    plot1 = figure.add_subplot(121)
    plot2 = figure.add_subplot(122)
    plot1.set_xlim(-4,4)
    plot1.set_ylim(-4,4)
    plot2.set_xlim(-4,4)
    plot2.set_ylim(-4,4)
    x0, y0, c= queue.get()
    scatter = plot1.scatter(
        x0,
        y0,
        c = c
    )    
    
    animation = animation.FuncAnimation(figure, update_plot, fargs = [queue,[scatter]], frames=100, interval = 50)
    pyplot.show()
    stop_event.set()
