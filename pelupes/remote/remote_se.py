from remote_worker import RemoteCodeInterface

from amuse.datamodel.parameters import ParametersMemento
from amuse.datamodel import Particles
from amuse.units import units

class forwarding_class_server(object):
    def __init__(self,_class, *arg,**kwarg):
        self.code=_class(*arg,**kwarg)
    def receive_parameters(self,parameters):
        self.code.parameters.reset_from_memento(parameters)
    def send_parameters(self):
        return self.code.parameters.copy()     
    def receive_particles(self,particles):        
        particles.synchronize_to(self.code.particles)
        channel=particles.new_channel_to(self.code.particles)
        channel.copy_all_attributes()
    def send_particles(self):
        return self.code.particles.copy()
    def evolve_model(self,tend):
        self.code.evolve_model(tend)
    def model_time(self):
        return self.code.model_time
        
class forwarding_class_client(object):
    def __init__(self, _class, server_arg=(),server_kwarg={},
                   code_arg=(),code_kwarg={}):

        try:
            self.remote=RemoteCodeInterface(*server_arg,**server_kwarg) 
        except Exception as ex:
            print "forwarding_class_client: startup of remote worker failed"
            raise ex

        self.start_remote(_class,*code_arg,**code_kwarg)

        self._parameters=self.request_parameters()
        self._particles=Particles()

        self.particles_accessed=True
        self.parameters_accessed=True
              
    @property  
    def particles(self):
        if not self.particles_accessed:
            self._particles=self.request_particles()
            self.particles_accessed=True
        return self._particles  

    @property
    def parameters(self):
        if not self.parameters_accessed:
            self._parameters=self.request_parameters()
            self.parameters_accessed=True
        return self._parameters

    def commit_particles(self):
        self.dispatch_particles(self._particles)      
        self.particles_accessed=False

    def commit_parameters(self):
        self.dispatch_particles(self._parameters)      
        self.particles_accessed=False

    def evolve_model(self,tend):
        if self.parameters_accessed:
          self.dispatch_parameters(self._parameters)  
        if self.particles_accessed:
          self.dispatch_particles(self._particles)
        self.remote_evolve_model(tend)
        self.particles_accessed=False

    def start_remote(self,_class, *arg, **kwarg):
        from remote_se import forwarding_class_server
        self.remote.assign("_class",_class)
        self.remote.assign("_arg",arg)
        self.remote.assign("_kwarg",kwarg)
        self.remote.assign("_server",forwarding_class_server)
        self.remote.execute("server=_server(_class,*_arg,**_kwarg)")

    def dispatch_parameters(self,parameters):
        self.remote.assign("_p", parameters)
        self.remote.execute("server.receive_parameters(_p)")
            
    def request_parameters(self):
        return self.remote.evaluate("server.send_parameters()")

    def dispatch_particles(self,particles):
        self.remote.assign("_p", particles)
        self.remote.execute("server.receive_particles(_p)")
            
    def request_particles(self):
        return self.remote.evaluate("server.send_particles()")

    def remote_evolve_model(self,tend):
        self.remote.assign("_t",tend)
        self.remote.execute("server.evolve_model(_t)")

    @property
    def model_time(self):
        return self.remote.evaluate("server.model_time()")

if __name__=="__main__":
    from amuse.community.sse.interface import SSE

    se=forwarding_class_client(SSE,server_kwarg=dict(redirection="none"))
    
    p=Particles(1, mass=1| units.MSun)
    
    print se.parameters
    print se.model_time
    
    se.particles.add_particles(p)
    se.evolve_model(1.| units.Gyr)
    print se.particles
    se.evolve_model(2.| units.Gyr)
    print se.particles
