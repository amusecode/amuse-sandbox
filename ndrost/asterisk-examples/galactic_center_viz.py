"""
Evolves a cluster in the potention of the galactic center
  
Uses the bridge integrator to couple different codes.
In this example a cluster is evolved  circling the galactic center, represented by a static potential.
"""

import numpy

from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system

from amuse.ext.bridge import bridge
from rotating_bridge import RotatingBridge,rotating_to_inertial,inertial_to_rotating

from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.fi.interface import Fi
from amuse.community.huayno.interface import Huayno
from amuse.community.gadget2.interface import Gadget2

from matplotlib import pyplot

from amuse.ic.kingmodel import new_king_model

from amuse.community.asterisk.interface import Asterisk

numpy.random.seed(1234567)

class GalacticCenterGravityCode(object):
    """
    Implements a code simulating the galactic center. As the center itself does
    not evolve we only need to define the 'get_gravity_at_point'
    and 'get_potential_at_point'. Note that both functions get arrays
    of points.
    """
    def __init__(self,R=1000.| units.parsec, M=1.6e10 | units.MSun, alpha=1.2):
        self.R=R
        self.M=M
        self.alpha=alpha

    def get_gravity_at_point(self,eps,x,y,z):
        r2=x**2+y**2+z**2
        r=r2**0.5
        m=self.M*(r/self.R)**self.alpha  
        fr=constants.G*m/r2
        ax=-fr*x/r
        ay=-fr*y/r
        az=-fr*z/r
        return ax,ay,az

    def get_potential_at_point(self,eps,x,y,z):
        r2=x**2+y**2+z**2
        r=r2**0.5
        c=constants.G*self.M/self.R**self.alpha    
        phi=c/(self.alpha-1)*(r**(self.alpha-1)-self.R**(self.alpha-1))
        return phi    

    def vcirc(self,r):  
        m=self.M*(r/self.R)**self.alpha  
        vc=(constants.G*m/r)**0.5
        return vc
        
def king_model_cluster(interface,N=1024,W0=3, Mcluster=4.e4 | units.MSun,
                                 Rcluster= .7 | units.parsec,parameters=[]):
    """
    helper function to setup an nbody king model cluster (returns code with particles)
    """      
    converter=nbody_system.nbody_to_si(Mcluster,Rcluster)

    parts=new_king_model(N,W0,convert_nbody=converter)
    parts.radius=0.02| units.parsec

    nb=interface(converter)
    for name,value in parameters:
      setattr(nb.parameters, name, value)

    nb.particles.add_particles(parts)

    return nb,converter

def shift_sys(system,dx,dy,dz,dvx,dvy,dvz,omega):
    """
    helper function to shift system
    """
    parts=system.particles.copy()
    parts.x=parts.x+dx
    parts.y=parts.y+dy
    parts.z=parts.z+dz
    parts.vx=parts.vx+dvx
    parts.vy=parts.vy+dvy
    parts.vz=parts.vz+dvz
#    parts=inertial_to_rotating(0.| units.Myr,omega,parts)
    
    channel=parts.new_channel_to(system.particles)
    channel.copy_attributes(["x","y","z","vx","vy","vz"])
    #      parts.copy_values_of_state_attributes_to(system.particles)
    
if __name__ in ('__main__', '__plot__'):

# parameter setup:
    N=1024
    W0=3
    Rinit=50. | units.parsec
    Mcluster=4.e4 | units.MSun
    Rcluster=1. | units.parsec
    nstep=500
    

# make cluster and Galactic center
    cluster,converter=king_model_cluster(Huayno,N,W0, Mcluster,Rcluster, parameters=[
                   ("epsilon_squared", (0.01 | units.parsec)**2), 
                                     # ("timestep",timestep/4)
                                      ] )
    center=GalacticCenterGravityCode()
    
    
# shift the cluster to an orbit around GC
# (note the systems share the same coordinate frame, although units may differ)    
    vcirc=center.vcirc(Rinit)

    omega=vcirc/Rinit

    shift_sys(cluster,Rinit,0| units.parsec,0.|units.parsec,
                  0| units.kms,0.*vcirc,0| units.kms,omega)

    Tperiod=2*numpy.pi/omega
    timestep=Tperiod/nstep
    
    print (Rcluster**3/constants.G/Mcluster)**0.5/timestep


# setup bridge; cluster is evolved under influence of GC
    sys=RotatingBridge(omega,verbose=False,use_threading=False)
    sys.add_system(cluster, (center,), False)   

    sys.time=0.| units.Myr

    time=[0.]
    de=[0.]
    E0=sys.kinetic_energy+sys.potential_energy+sys.jacobi_potential_energy

# evolve and make plots
    times=Tperiod*numpy.array(range(nstep+1))/(1.*nstep)

    viz=Asterisk(converter,redirection="none",channel_type="sockets")
    viz.initialize_code()
    
    tmp=cluster.particles.copy()
    tmp.red=1.0
    tmp.green=0.55
    tmp.blue=0.35
    tmp.x=tmp.x-Rinit
    
    viz.particles.add_particles(tmp)
    viz.store_view("0 Myr")

    for i,t in enumerate(times):
        sys.evolve_model(t,timestep=timestep)

        E=sys.kinetic_energy+sys.potential_energy+sys.jacobi_potential_energy
    
        time.append(t.number)

        de.append(abs((E0-E)/E0))
        tnow=sys.model_time
        print tnow

        tmp=sys.particles.copy()
        tmp.x=tmp.x-Rinit
        channel=tmp.new_channel_to(viz.particles)
        channel.copy_attributes(["x","y","z"])
        viz.store_view(t.in_(units.Myr))
        

    viz.stop()
     
    print max(de)        
    
