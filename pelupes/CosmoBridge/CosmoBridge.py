from amuse.ext.bridge import bridge

from amuse.units import quantities

class VacuumBoundary_CosmoBridge(bridge):
    def __init__(self,z0,cosmology,*args,**kwargs):
        self.cosmology=cosmology
        bridge.__init__(self,*args,**kwargs)
        self.time=self.cosmology.agefromz(z0)
        if 'approx' in kwargs:
          self.approx=kwargs['approx']
        else:
          self.approx=True

    def cosmo_fac(self,time):
        a=self.cosmology.afromage(time)
        d2adt=self.cosmology.d2adtau(a)*self.cosmology.hubble0**2        
        acc_omegam=0.5*self.cosmology.hubble0**2*self.cosmology.omegam/a**3
        return d2adt+ acc_omegam    

    def cosmo_fac_approx(self,time):
        return self.cosmology.omegal*self.cosmology.hubble0**2
        
    def cosmo_kick(self,system,partners,dt,cosmo_acc_fac):
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
        
        ax+=parts.x*cosmo_acc_fac
        ay+=parts.y*cosmo_acc_fac
        az+=parts.z*cosmo_acc_fac
        
        parts.vx=parts.vx+ax*dt
        parts.vy=parts.vy+ay*dt
        parts.vz=parts.vz+az*dt

        channel=parts.new_channel_to(system.particles)
        channel.copy_attributes(["vx","vy","vz"])                          
        if(self.verbose):
            print ".. done"

    def kick_systems(self,dt):

        if self.approx:
            cosmo_acc_fac=self.cosmo_fac_approx(self.time)
        else:
            cosmo_acc_fac=(self.cosmo_fac(self.time)+self.cosmo_fac(self.time+dt))/2

        for x in self.systems:
            if self.do_sync[x]:
                if hasattr(x,"synchronize_model"):
                    if(self.verbose): print x.__class__.__name__,"is synchronizing",
                    x.synchronize_model()    
                    if(self.verbose):  print ".. done"
        for x in self.systems:
            if hasattr(x,"particles"):
                self.cosmo_kick(x,self.partners[x],dt,cosmo_acc_fac)                  

        self._kick_time+=dt
        return 0
    
    @property 
    def cosmo_energy_approx(self):
        r2=self.particles.lengths()
        return (-0.5*self.cosmology.omegal*self.cosmology.hubble0**2*r2).sum()


class PeriodicBoundary_CosmoBridge(bridge):
    def __init__(self,z0,cosmology,*args,**kwargs):
        self.cosmology=cosmology
        bridge.__init__(self,*args,**kwargs)
        self.time=self.cosmology.agefromz(z0)

    def cosmo_fac(self,time):
        a=self.cosmology.afromage(time)
        d2adt=self.cosmology.d2adtau(a)*self.cosmology.hubble0**2
        return d2adt   

    def cosmo_kick(self,system,partners,dt,cosmo_acc_fac):
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
        
        ax+=parts.x*cosmo_acc_fac
        ay+=parts.y*cosmo_acc_fac
        az+=parts.z*cosmo_acc_fac
        
        parts.vx=parts.vx+ax*dt
        parts.vy=parts.vy+ay*dt
        parts.vz=parts.vz+az*dt

        channel=parts.new_channel_to(system.particles)
        channel.copy_attributes(["vx","vy","vz"])                          
        if(self.verbose):
            print ".. done"

    def kick_systems(self,dt):

        cosmo_acc_fac=(self.cosmo_fac(self.time)+self.cosmo_fac(self.time+dt))/2

        for x in self.systems:
            if self.do_sync[x]:
                if hasattr(x,"synchronize_model"):
                    if(self.verbose): print x.__class__.__name__,"is synchronizing",
                    x.synchronize_model()    
                    if(self.verbose):  print ".. done"
        for x in self.systems:
            if hasattr(x,"particles"):
                self.cosmo_kick(x,self.partners[x],dt,cosmo_acc_fac)                  

        self._kick_time+=dt
        return 0
    
