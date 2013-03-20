"""
   convention: positive omega means: rotating frame is rotating counterclockwise
  (i.e. in the rotating frame a particle at rest in the inertial frame seems to
  rotate clockwise)

  non-canonical coordinates
  
  get_gravity_at_point should return accelerations in rotating frame

"""

from amuse.ext.bridge import bridge

from amuse.units import quantities

from numpy import cos,sin

def inertial_to_rotating(t,omega,parts):
  x=parts.x
  y=parts.y
  vx=parts.vx
  vy=parts.vy
  rotating=parts.copy()
  rotating.x=x*cos(omega*t)+y*sin(omega*t)
  rotating.y=-x*sin(omega*t)+y*cos(omega*t)
  rotating.vx=(vx+y*omega)*cos(omega*t)+(vy-x*omega)*sin(omega*t)
  rotating.vy=-(vx+y*omega)*sin(omega*t)+(vy-x*omega)*cos(omega*t)
  return rotating
  
def rotating_to_inertial(t,omega,parts):     
  return inertial_to_rotating(t,-omega,parts)

class RotatingBridge(bridge):
    def __init__(self,omega,*args,**kwargs):
        self.omega=omega
        bridge.__init__(self,*args,**kwargs)

    def rotating_kick(self,system,partners,dt):
        parts=system.particles.copy()
        ax=quantities.zero
        ay=quantities.zero
        az=quantities.zero
        if(self.verbose):  
            print system.__class__.__name__,"receives kick from",
        for y in partners:
            if system is not y:
                if(self.verbose):  
                    print y.__class__.__name__,
                _ax,_ay,_az=y.get_gravity_at_point(parts.radius,parts.x,parts.y,parts.z)
                ax+=_ax
                ay+=_ay
                az+=_az
        if self.omega.value_in(self.omega.unit) != 0.:
          vx0=parts.vx
          vy0=parts.vy
          omega=2*self.omega
          a1_omega=(ax+self.omega**2*parts.x)/omega
          a2_omega=(ay+self.omega**2*parts.y)/omega
          parts.vx=(vx0-a2_omega)*cos(omega*dt)+(vy0+a1_omega)*sin(omega*dt)+a2_omega
          parts.vy=-(vx0-a2_omega)*sin(omega*dt)+(vy0+a1_omega)*cos(omega*dt)-a1_omega
          parts.vz=parts.vz+az*dt
        else:
          parts.vx=parts.vx+ax*dt
          parts.vy=parts.vy+ay*dt
          parts.vz=parts.vz+az*dt
        channel=parts.new_channel_to(system.particles)
        channel.copy_attributes(["vx","vy","vz"])                          
        if(self.verbose):
            print ".. done"

    def kick_systems(self,dt):
        for x in self.systems:
            if self.do_sync[x]:
                if hasattr(x,"synchronize_model"):
                    if(self.verbose): print x.__class__.__name__,"is synchronizing",
                    x.synchronize_model()    
                    if(self.verbose):  print ".. done"
        for x in self.systems:
            if hasattr(x,"particles"):
                self.rotating_kick(x,self.partners[x],dt)                  
        return 0
 
    def get_effective_potential_at_point(self,radius,x,y,z):
        err=0
        pot=bridge.get_potential_at_point(self,radius,x,y,z)
        r2=x**2+y**2+z**2
        pot+=-0.5*self.omega**2*r2
        return pot
        
    @property
    def jacobi_potential_energy(self):
        parts=self.particles
        return -0.5*(parts.mass*self.omega**2*(parts.x**2+parts.y**2)).sum()
       
       
