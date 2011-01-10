import numpy
import scipy.interpolate
import cPickle

class tidal_field(object):
  def __init__(self, tensor):
    self.tidal_tensor=numpy.array(tensor,'d')
      
  def get_gravity_at_point(self,eps,x,y,z):    
    pos=numpy.array( (0.,) * 3, 'd' )
    acc=numpy.array( (0.,) * 3, 'd' )    
    ax=numpy.array(eps,'d')
    ay=numpy.array(eps,'d')
    az=numpy.array(eps,'d')
    pos[0]=x
    pos[1]=y
    pos[2]=z
    acc=numpy.dot(self.tidal_tensor,pos)
    ax.fill(acc[0])
    ay.fill(acc[1])
    az.fill(acc[2])
    return ax,ay,az,0

  def get_potential_at_point(self,eps,x,y,z):    
    pos=numpy.array( (0.,) * 3, 'd' )
    phi=numpy.array(eps,'d')
    pos[0]=x
    pos[1]=y
    pos[2]=z
    phi.fill(-0.5*numpy.dot(pos,numpy.dot(self.tidal_tensor,pos)))
    return phi,0

class time_dependent_tidal_field(tidal_field):
  def __init__(self, time, tides):
    print "tmin,tmax:", time[0],time[-1] 
    self.f=scipy.interpolate.interp1d(time-time[0],tides,axis=0)
    self.tidal_tensor=self.f(time[0])
    
  def evolve(self,t):
    self.tidal_tensor=self.f(t)  
      


