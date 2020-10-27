import cPickle

from intercode import Code

def test():
    instance=Code(channel_type="mpi",redirection="none")

    port=instance.internal__open_port()

    print port
    
    instance.stop()

def test2():
    instance1=Code(channel_type="sockets",redirection="none",number_of_workers=4)
    instance2=Code(channel_type="sockets",redirection="none")

    port=instance1.internal__open_port()

    print port

    id1 = instance1.internal__accept_on_port.asynchronous(port)
    id2 = instance2.internal__connect_to_port.asynchronous(port)

    print id1.result(),id2.result()

    print "here"

    instance1.stop()
    instance2.stop()

def test3():
    number_of_workers=2
    instance1=Code(channel_type="sockets",redirection="none",number_of_workers=number_of_workers)
    instance2=Code(channel_type="mpi",redirection="none",number_of_workers=number_of_workers)

    port=instance1.internal__open_port()

    id1 = instance1.internal__accept_on_port.asynchronous(port)
    id2 = instance2.internal__connect_to_port.asynchronous(port)

    print instance1.grid

    instance1.grid.dens=instance1.grid.x**2

    print instance1.grid.dens

    instance1.send_.asynchronous(id1.result())
    instance2.receive_.asynchronous(id2.result())

    print instance2.grid.dens

    instance1.stop()
    instance2.stop()


if __name__=="__main__":
    test3()
