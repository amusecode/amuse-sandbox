""" Implementation of gas cooling functions for a set of gas particles. """
import numpy

from amuse.units import units, quantities, constants

class Cooling(object):
    def __init__(self, particles):
        self.particles = particles

        self.cooling_function = self.my_cooling_function
        self.heating_function = self.gerritsen_heating_function

    def evolve_model(self, model_time):
        self.particles.u = self.evolve_internal_energy(self.particles.u, model_time, self.particles.rho/self.mu())

    def evolve_internal_energy(self, u_0, dt, n_H, du_dt_adiabatic=quantities.zero):
        def full_function(u):
            return du_dt_adiabatic * u/u_0 + (self.heating_function() - n_H * self.cooling_function(self.T_from_u(u))) / self.mu()

        u_out = self.integrate_ode(full_function, u_0, dt)

        return u_out

    def integrate_ode(self, function, x, t_end, eps=0.01):
        """
        Integrates the given ordinary differential equation of the form:
        dx/dt = function(x)
        for a time 't_end', using the initial value 'x'.
        The routine takes small steps, such that (abs(dx) <= eps * x)
        """
        t = 0 * t_end
        while t < t_end:
            fx = function(x)
            step = min( (t_end-t), ((eps*x)/abs(fx)).amin() )
            t += step
            x += fx * step
        return x

    def T_from_u(self, u):
        return 2.0/3.0 * u * self.mu() / constants.kB

    def mu(self, X=None, Y=0.25, Z=0.02, x_ion=0.1):
        """
            Compute the mean molecular weight in kg (the average weight of particles in a gas)
            X, Y, and Z are the mass fractions of Hydrogen, of Helium, and of metals, respectively.
            x_ion is the ionisation fraction (0 < x_ion < 1), 1 means fully ionised
        """
        if X is None:
            X = 1.0 - Y - Z
        elif abs(X + Y + Z - 1.0) > 1e-6:
            raise Exception("Error in calculating mu: mass fractions do not sum to 1.0")
        return constants.proton_mass / (X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0 + Z*x_ion/2.0)

    def gerritsen_heating_function(self, G_0=10, eps=0.05):
        return 10.0**-24 * eps * G_0 | units.erg / units.s

    def no_heating_function(self):
        return quantities.zero
    def gerritsen_cooling_function(self, T, logT= None, a=3.24, b=0.170): # x=1e-1
        if logT is None:
            logT = numpy.log10(T.value_in(units.K))
        condlist = [logT <= 6.2, logT >= 6.2]
        choicelist = [10.0**-21.0 * (10**(-0.1-1.88*(5.23-logT)**4) + 10**(-a-b*(4-logT)**2)), 10.0**-22.7]
        return (units.erg*units.cm**3/units.s).new_quantity(numpy.select(condlist, choicelist))

    def my_cooling_function(self, T, logT=None, a=3.24, b=0.170): # x=1e-1
        if logT is None:
            logT = numpy.log10(T.value_in(units.K))
        condlist = [logT <= 6.2, logT >= 6.2]
        choicelist = [10.0**-21.0 * (10**(-0.1-1.88*(5.23-logT)**4) + 10**(-a-b*abs(4-logT)**3)), 10.0**-22.7]
        return (units.erg*units.cm**3/units.s).new_quantity(numpy.select(condlist, choicelist))
