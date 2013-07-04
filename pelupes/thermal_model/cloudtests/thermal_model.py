import numpy

from matplotlib import pyplot

from amuse.units import units,constants
from amuse.units.quantities import VectorQuantity,zero
from amuse.datamodel import ParticlesOverlay

class SimplifiedThermalModel(object):
  def __init__(self, n0=0.05 | units.cm**-3, T0=1.e4 | units.K,
                   Tmin=20 | units.K, alpha=5., reference_heating=1.e-25 | units.erg/units.s):
      self.reference_mu=(2.2 | units.amu)
      self.rho0=n0*self.reference_mu
      self.T0=T0
      self.Tmin=Tmin
      self.alpha=alpha
      self.reference_heating=reference_heating

  def equilibrium_temperature(self,rho):
    xclip=(rho/self.rho0)
    return self.Tmin+ (self.T0-self.Tmin)/(1.+numpy.log10(1.+xclip)**self.alpha)

  def mu(self,rho):
    return numpy.ones_like(rho.number)*self.reference_mu

  def gamma(self,rho):
    return numpy.ones_like(rho.number)*(self.reference_heating)

  def equilibrium_u(self,rho):
    return constants.kB*self.equilibrium_temperature(rho)/self.mu(rho)
    
  def tau(self,rho):
    return (constants.kB*self.equilibrium_temperature(rho)/self.gamma(rho))

  def evolve_u(self,dt,rho,u0,dudt=None):
    u_eq=self.equilibrium_u(rho)
    tau=self.tau(rho)
    if dudt is not None:
      condition1= 1.*(dudt*tau < (u0-u_eq))
      condition2= 1.-condition1
      fac=1./numpy.maximum(1-dudt/u0*tau ,1.e-5)
      u_eq=(u_eq*fac)*condition1+(u_eq+dudt*tau)*condition2
      tau=(tau*fac)*condition1+tau*condition2
    return u_eq+(u0-u_eq)*numpy.exp(-dt/tau)

  def evolve_u_radiated_energy(self,dt,rho,u0,dudt=None):
    u_eq=self.equilibrium_u(rho)
    tau=self.tau(rho)
    if dudt is not None:
      condition1= 1.*(dudt*tau < (u0-u_eq))
      condition2= 1.-condition1
      fac=1./numpy.maximum(1-dudt/u0*tau ,1.e-5)
      u_eq=(u_eq*fac)*condition1+(u_eq+dudt*tau)*condition2
      tau = (tau*fac)*condition1+tau*condition2
    u1=u_eq+(u0-u_eq)*numpy.exp(-dt/tau)
    deltau=u1-u0
    rad=-deltau
    if dudt is not None:
      condition_eq= 1.-(u0!=u_eq)*(u1!=u0)
      logfac=numpy.log( ( (u1-u_eq)+u_eq*condition_eq)/ ((u0-u_eq)+u_eq*condition_eq) )
      rad=(-deltau*fac - dudt*tau*(u_eq/u0)*logfac)*condition1 +  \
          (-deltau-dudt*tau*logfac)*condition2+ \
          (dudt*dt)*condition_eq
    return u1,rad
# rad>0 -> cooling 
# rad<0 -> heating

class SPH_with_Thermal_Model(object):
    def __init__(self,sph_code,thermal_model,timestep, calc_net_luminosity=True):
        self.sph_code=sph_code()
        self.thermal_model=thermal_model()
        self.timestep=timestep
        self.calc_net_luminosity=calc_net_luminosity
        self._gas_particles=ParticlesOverlay(self.sph_code.gas_particles)
        if self.calc_net_luminosity:
          self.radiated_energy=zero
          self.total_luminosity=zero
    
    @property
    def model_time(self):
      return self.sph_code.model_time
    
    @property
    def parameters(self):
        return self.sph_code.parameters
        
    @property
    def particles(self):
        return self.sph_code.particles    
    @property
    def star_particles(self):
        return self.sph_code.star_particles    
    @property
    def dm_particles(self):
        return self.sph_code.dm_particles
    @property
    def gas_particles(self):
        return self._gas_particles    
    @property
    def sink_particles(self):
        return self.sph_code.sink_particles    

    def evolve_u(self,dt):
        density=self.gas_particles.density
        u0=self.gas_particles.u
        if self.parameters.isothermal_flag:
          dudt=self.gas_particles.du_dt
        else:
          dudt=None
        if self.calc_net_luminosity:
          u,lum=self.thermal_model.evolve_u_radiated_energy(dt,density,u0,dudt=dudt)
          self.gas_particles.u=u
          self.gas_particles.specific_net_luminosity=lum/dt
          self.total_luminosity=(self.gas_particles.mass*lum).sum()/dt
          self.radiated_energy+=(self.gas_particles.mass*lum).sum()
        else:
          self.gas_particles.u=self.thermal_model.evolve_u(dt,density,u0,dudt=dudt)

    def evolve_model(self,tend):
        tnow=self.model_time
        while tnow < (tend - self.timestep/2):
          dt=tend-self.model_time
          if self.timestep<dt:
            dt=self.timestep
          print "1st evolve_u",
          self.evolve_u(dt/2)
          print "evolve sph",
          self.sph_code.evolve_model(tnow+dt)
          print "2nd evolve_u",
          self.evolve_u(dt/2)
          print
          tnow=self.model_time


def attach_thermal_model(baseclass):
    class newclass(baseclass):
      def __init__(self,*args, **kwargs):
          self.thermal_model=kwargs.pop("thermal_model")()
          self.timestep=kwargs.pop("timestep")
          self.calc_net_luminosity=kwargs.pop("calc_net_luminosity")
          baseclass.__init__(self, *args, **kwargs)
          self._gas_particles=ParticlesOverlay(self.overridden().gas_particles)
          if self.calc_net_luminosity:
              self.radiated_energy=zero
              self.total_luminosity=zero

      @property
      def gas_particles(self):
          return self._gas_particles    

      def evolve_u(self,dt):
          density=self.gas_particles.density
          u0=self.gas_particles.u
          if self.parameters.isothermal_flag:
            dudt=self.gas_particles.du_dt
          else:
            dudt=None
          if self.calc_net_luminosity:
            u,lum=self.thermal_model.evolve_u_radiated_energy(dt,density,u0,dudt=dudt)
            self.gas_particles.u=u
            self.gas_particles.specific_net_luminosity=lum/dt
            self.total_luminosity=(self.gas_particles.mass*lum).sum()/dt
            self.radiated_energy+=(self.gas_particles.mass*lum).sum()
          else:
            self.gas_particles.u=self.thermal_model.evolve_u(dt,density,u0,dudt=dudt)

      def evolve_model(self,tend):
          tnow=self.model_time
          while tnow < (tend - self.timestep/2):
            dt=tend-self.model_time
            if self.timestep<dt:
              dt=self.timestep
            print "1st evolve_u",
            self.evolve_u(dt/2)
            print "evolve sph",
            try:
              baseclass.evolve_model(self,tnow+dt)
            except AttributeError:
              self.overridden().evolve_model(tnow+dt)  
            print "2nd evolve_u",
            self.evolve_u(dt/2)
            print
            tnow=self.model_time
    return newclass

if __name__=="__main__":
    pass
  
    
