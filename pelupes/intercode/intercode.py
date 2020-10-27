import numpy
from mpi4py import MPI

from amuse.support.interface import InCodeComponentImplementation

from amuse.rfi.core import PythonCodeInterface, remote_function

from amuse.datamodel import CartesianGrid

class CodeImplementation(object):
    
    def __init__(self):
        self.comm=MPI.COMM_WORLD
        self.myrank=self.comm.Get_rank()
        self.N=self.comm.Get_size()
        self.Ngrid=3*4*5
        n=self.Ngrid/self.N
        x = (numpy.arange(n)+self.myrank*n)/(1.*self.Ngrid)
        self.local_imin=self.myrank*n
        self.local_imax=(self.myrank+1)*n-1
        self.dens = numpy.zeros(n)
        
    def get_range(self,imin,imax):
        imin.value=0
        imax.value=self.Ngrid-1
        return 0
        
    def get_x(self,index,x):
        x.value=index/(1.*self.Ngrid)
        return 0

    def get_dens(self,index,dens,N):
        a=(index>=self.local_imin)*(index<=self.local_imax)
        _dens=numpy.zeros(N)
        _dens[a]=self.dens[index[a]-self.local_imin]
        dens.value=numpy.zeros(N)
        _dens=self.comm.Reduce(_dens, dens.value, MPI.SUM,root=0)
        return 0

    def set_dens(self,index,dens,N):
        a=(index>=self.local_imin)*(index<=self.local_imax)
        self.dens[index[a]-self.local_imin]=dens[index[a]]
        return 0

    def send_(self, comm_id):
        comm=self._interface.communicators[comm_id]
        rank=comm.Get_rank()
        comm.Send(self.dens, dest=rank, tag=rank)
        return 0

    def receive_(self,comm_id):
        comm=self._interface.communicators[comm_id]
        rank=comm.Get_rank()
        comm.Recv(self.dens, source=rank, tag=rank)
        return 0


    def parallel_get_dens(self,index, dens, N):
        a=(index>=self.local_imin)*(index<=self.local_imax)
        _dens=numpy.zeros(N)
        _dens[a]=self.dens[index[a]-self.local_imin]
        dens.value=_dens
        return 0
    

class CodeInterface(PythonCodeInterface):
    
    def __init__(self, **options):
        PythonCodeInterface.__init__(self, implementation_factory = CodeImplementation, **options)

    @remote_function
    def get_range():
        returns(imin=0,imax=0)

    @remote_function
    def get_x(index=0):
        returns(x=0.)
        
    @remote_function(must_handle_array=True)
    def get_dens(index=0):
        returns(dens=0.)

    @remote_function(must_handle_array=True)
    def set_dens(index=0, dens=0.):
        returns()

    @remote_function
    def send_(comm_id=0):
        returns()

    @remote_function
    def receive_(comm_id=0):
        returns()

    
    #~ @parallel_function
    #~ def parallel_get_dens(index=0):
        #~ returns(dens=0.)

class Code(InCodeComponentImplementation):
    
    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self, CodeInterface(**options), **options)

    def define_grids(self, object):        
        object.define_grid('grid',axes_names = ['x'],grid_class=CartesianGrid)
        object.set_grid_range('grid', 'get_range')
        object.add_getter('grid', 'get_dens', names=('dens',))
        object.add_setter('grid', 'set_dens', names=('dens',))
        object.add_getter('grid', 'get_x', names=('x',))

