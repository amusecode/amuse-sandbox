import numpy

from matplotlib import pyplot

from amuse.units import units,constants
from amuse.units.quantities import VectorQuantity,zero
from amuse.datamodel import ParticlesOverlay

from thermal_model import SimplifiedThermalModel

def plot_eq_curves():
    n=10**(-2 + 5* numpy.arange(1001)/1000. ) | units.cm**-3
    

    f=pyplot.figure(figsize=(8,6))

    for n0 in [0.025 | units.cm**-3,0.05 | units.cm**-3,0.1 | units.cm**-3]:
      thm=SimplifiedThermalModel(n0=n0)
      rho=n*thm.reference_mu
      t=thm.equilibrium_temperature(rho)
      pyplot.loglog(n.number,t.number,lw=2.)    
    pyplot.ylim(10,1.e5)
    pyplot.ylabel("T (K)")
    pyplot.xlabel("n ("+str(n.unit)+")")
    pyplot.savefig("n_vs_T_eq.eps")
    pyplot.show()

    f=pyplot.figure(figsize=(8,6))
    for n0 in [0.025 | units.cm**-3,0.05 | units.cm**-3,0.1 | units.cm**-3]:
      thm=SimplifiedThermalModel(n0=n0)
      rho=n*thm.reference_mu
      t=thm.equilibrium_temperature(rho)    
      pyplot.loglog(n.number,(n*t).value_in(units.K*units.cm**-3),lw=2.)    
    pyplot.xlabel("n ("+str(n.unit)+")")
    pyplot.ylabel("P ("+str(units.K*units.cm**-3)+")")
    pyplot.savefig("n_vs_P_eq.eps")
    pyplot.show()

    f=pyplot.figure(figsize=(8,6))
    for n0 in [0.025 | units.cm**-3,0.05 | units.cm**-3,0.1 | units.cm**-3]:
      thm=SimplifiedThermalModel(n0=n0)
      rho=n*thm.reference_mu
      t=thm.equilibrium_temperature(rho)
      tau=thm.tau(rho).in_(units.Myr)
      pyplot.loglog(n.number,tau.number,lw=2.)
    pyplot.xlabel("n ("+str(n.unit)+")")
    pyplot.ylabel("tau ("+str(tau.unit)+")")
    pyplot.savefig("n_vs_tau.eps")
    pyplot.show()
  
def plot_time_evolution_lum():  
    thm=SimplifiedThermalModel()

    rho=([100] | units.cm**-3)*thm.reference_mu
    
    u_eq=thm.equilibrium_u(rho)
    
    f=pyplot.figure(figsize=(8,6))

    tau=thm.tau(rho)
    refdudt=u_eq/tau

    c="rgbyc"

    dt=.1*thm.tau(rho)
    
    for i,dudt in enumerate([-refdudt/2.,0.*refdudt,None,refdudt,2*refdudt]):    
        
#        u=u_eq+20*refdudt*tau
        u=10.*u_eq
    
        tnow=dt*0.
        time=[tnow]
        ltime=[tnow]
        temp=[thm.reference_mu*u/constants.kB]
        if dudt is None:
          lum=[(u-u_eq)/tau]
          teq=thm.reference_mu*u_eq/constants.kB
        else:
          lum=[(u-u_eq)/tau]
          teq=thm.reference_mu*(u_eq+dudt*tau)/constants.kB
          
        while tnow< 10*thm.tau(rho):
          u,l=thm.evolve_u_radiated_energy(dt,rho,u,dudt=dudt)
          tnow=tnow+dt
          time.append(tnow)
          temp.append(thm.reference_mu*u/constants.kB)
          ltime.append(tnow-dt/2)
          lum.append(l/dt)
                  
        time=VectorQuantity.new_from_scalar_quantities(*time)
        ltime=VectorQuantity.new_from_scalar_quantities(*ltime)
        temp=VectorQuantity.new_from_scalar_quantities(*temp)  
        lum=VectorQuantity.new_from_scalar_quantities(*lum)
                
        pyplot.subplot(211)      
        pyplot.semilogy(time/tau,temp.number,c[i],lw=2.)
        pyplot.semilogy(time/tau,numpy.ones_like(time.number)*teq.number,c[i]+':',lw=2.)
        pyplot.xlabel("time/tau")
        pyplot.ylabel("T (K)")
        
        pyplot.ylim(1,1000)
        pyplot.xlim(0,10)
        pyplot.subplot(212)      
        lum=lum.in_(units.milli(units.W)/units.kg)
        pyplot.plot(ltime/tau,lum.number,c[i],lw=2.)
        if dudt is not None:
          pyplot.plot(time/tau,(numpy.ones_like(time.number)*dudt).value_in(lum.unit),c[i]+':',lw=2.)
        pyplot.xlabel("time/tau")
        pyplot.ylabel("L ("+str(lum.unit)+")")
        pyplot.xlim(0,10)
#        pyplot.ylim(-1.e-18,1.e-18)
    
    pyplot.savefig("T_vs_time_wL.eps")
#    pyplot.show()
  
def plot_time_evolution():
    thm=SimplifiedThermalModel()
    
    rho=([100] | units.cm**-3)*thm.reference_mu
    
    u_eq=thm.equilibrium_u(rho)
    tau=thm.tau(rho)
    refdudt=u_eq/tau
    
    c="rgbycm"

    f=pyplot.figure(figsize=(8,6))
    
    i=-1
    for dt in [2.*tau,0.5*tau,0.05*tau]:    
      for u in [0.2*u_eq,10*u_eq]:

          i=(i+1)%len(c)
          tnow=dt*0.
          time=[tnow]
          temp=[thm.reference_mu*u/constants.kB]
          
          teq=thm.reference_mu*(thm.equilibrium_u(rho))/constants.kB
          
          while tnow< 10*tau-dt/2:
#            u,l=thm.evolve_u_radiated_energy(dt,rho,u)
            u=thm.evolve_u(dt,rho,u)
            tnow=tnow+dt
            time.append(tnow)
            temp.append(thm.reference_mu*u/constants.kB)
            
          time=VectorQuantity.new_from_scalar_quantities(*time)
          temp=VectorQuantity.new_from_scalar_quantities(*temp)  
          
          pyplot.semilogy(time/tau,temp.number,c[i],lw=2.)
    pyplot.semilogy(time/tau,numpy.ones_like(time.number)*teq.number,c[i]+':',lw=2.)

    pyplot.xlabel("time/tau")
    pyplot.ylabel("T (K)")

    pyplot.savefig("t_vs_Teq.eps")
 
    pyplot.show()

  
  
if __name__=="__main__":
    
    pyplot.ion()

    plot_eq_curves()    
    plot_time_evolution()
    plot_time_evolution_lum()
    
    raw_input()
