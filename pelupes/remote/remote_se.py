from remote_worker import RemoteCodeInterface

from amuse.datamodel.parameters import ParametersMemento
from amuse.datamodel import Particles
from amuse.units import units

class forwarding_class_server(object):
    def __init__(self,_class):
      self.code=_class(channel_type="sockets")
    def receive_parameters(self,parameters):
        self.code.parameters.reset_from_memento(parameters)
    def send_parameters(self):
        return self.code.parameters.copy()     
    def receive_particles(self,particles):
        add_set=particles.difference(self.code.particles)
        remove_set=self.code.particles.difference(particles)

        if len(remove_set)>0: self.code.particles.remove_particles(remove_set)
        if len(add_set)>0: self.code.particles.add_particles(add_set)
        self.code.commit_particles() 

        channel=particles.new_channel_to(self.code.particles)
        channel.copy_all_attributes()
    def send_particles(self):
        return self.code.particles.copy()
    def evolve_model(self,tend):
        self.code.evolve_model(tend)
        
class forwarding_class_client(object):
    def __init__(self, _class, *arg, **kwargs):

        try:
          self.remote=RemoteCodeInterface(*arg,**kwargs) 
        except Exception as ex:
          print "forwarding_class_client: startup of remote worker failed"
          raise ex

        self.start_remote(_class)

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

    def start_remote(self,_class):
        from remote_se import forwarding_class_server
        self.remote.assign("_class",_class)
        self.remote.assign("_server",forwarding_class_server)
        self.remote.execute("server=_server(_class)")

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

if __name__=="__main__":
    from amuse.community.sse.interface import SSE

    se=forwarding_class_client(SSE,redirection="none")
    
    p=Particles(1, mass=1| units.MSun)
    
    print se.parameters
    
    se.particles.add_particles(p)
    se.evolve_model(1.| units.Gyr)
    print se.particles
    se.evolve_model(2.| units.Gyr)
    print se.particles
