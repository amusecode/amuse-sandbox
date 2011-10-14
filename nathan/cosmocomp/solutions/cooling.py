import numpy

try:
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
    from amuse.plot import plot, xlabel, ylabel, loglog
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.console import set_printing_strategy
from amuse.units import units
from amuse.units import constants
from amuse.units.quantities import zero


def gerritsen_heating_function(G_0 = 10, eps = 0.05):
    return 10.0**-24 * eps * G_0 | units.erg / units.s

def gerritsen_cooling_function(T, logT = None, a = 3.24, b = 0.170): # x=1e-1
    if logT is None:
        logT = numpy.log10(T.value_in(units.K))
    condlist = [logT <= 6.2, logT >= 6.2]
    choicelist = [10.0**-21.0 * (10**(-0.1-1.88*(5.23-logT)**4) + 10**(-a-b*(4-logT)**2)), 10.0**-22.7]
    return (units.erg*units.cm**3/units.s).new_quantity(numpy.select(condlist, choicelist))

def my_cooling_function(T, logT = None, a = 3.24, b = 0.170): # x=1e-1
    if logT is None:
        logT = numpy.log10(T.value_in(units.K))
    condlist = [logT <= 6.2, logT >= 6.2]
    choicelist = [10.0**-21.0 * (10**(-0.1-1.88*(5.23-logT)**4) + 10**(-a-b*abs(4-logT)**3)), 10.0**-22.7]
    return (units.erg*units.cm**3/units.s).new_quantity(numpy.select(condlist, choicelist))


def plot_cooling_function(cooling_function, figure_name):
    if not HAS_MATPLOTLIB:
        return
    points_per_decade = 20
    logT = numpy.linspace(1.0, 8.0, num = points_per_decade * 7 + 1)
    T = units.K.new_quantity(10**logT)
    
    pyplot.figure(figsize = (12, 10))
    pyplot.gca().set_ylim([1.0e-28, 1.0e-21])
    loglog(T, cooling_function(T))
    xlabel('log(T)')
    ylabel('log($\Lambda$)')
    pyplot.savefig(figure_name)
    print "\nPlot of cooling function was saved to: ", figure_name
    pyplot.close()

def plot_cooling_vs_heating_n(n_H, cooling_function, heating_function, figure_name):
    if not HAS_MATPLOTLIB:
        return
    points_per_decade = 20
    logT = numpy.linspace(1.0, 8.0, num = points_per_decade * 7 + 1)
    T = units.K.new_quantity(10**logT)
    cool = cooling_function(T)
    heat = heating_function() * numpy.ones_like(logT)
    
    pyplot.figure(figsize = (12, 10))
    loglog(T, heat, label = "$\Gamma$")
    loglog(T, n_H * cool, label = "$n_H \Lambda$")
    xlabel('log(T)')
    ylabel('log($n_H \Lambda$), log($\Gamma$)')
    pyplot.legend(loc=4)
    pyplot.savefig(figure_name)
    print "\nPlot of cooling/heating at constant density was saved to: ", figure_name
    pyplot.close()

def plot_cooling_vs_heating_P(P, cooling_function, heating_function, figure_name):
    if not HAS_MATPLOTLIB:
        return
    points_per_decade = 20
    logT = numpy.linspace(1.0, 8.0, num = points_per_decade * 7 + 1)
    T = units.K.new_quantity(10**logT)
    n_H = P / T
    cool = cooling_function(T)
    heat = heating_function() * numpy.ones_like(logT)
    
    pyplot.figure(figsize = (12, 10))
    loglog(T, heat, label = "$\Gamma$")
    loglog(T, n_H * cool, label = "$n_H \Lambda$")
    xlabel('log(T)')
    ylabel('log($n_H \Lambda$), log($\Gamma$)')
    pyplot.legend(loc=4)
    pyplot.savefig(figure_name)
    print "\nPlot of cooling/heating at constant pressure was saved to: ", figure_name
    pyplot.close()


def plot_pressure_density_in_thermal_equilibrium(cooling_function, heating_function, figure_name):
    if not HAS_MATPLOTLIB:
        return
    points_per_decade = 100
    log_n_H = numpy.linspace(-2.0, 3.0, num = points_per_decade * 5 + 1)
    n_H = (units.cm**-3).new_quantity(10**log_n_H)
    
    logT_equi = numpy.array([log_equilibrium_temperature(dens, cooling_function, 
        heating_function, logT_min = -1.0, logT_max = 5.0) for dens in n_H])
    T_equi = 10.0**logT_equi | units.K
    P = (n_H * T_equi).as_quantity_in(units.K * units.cm**-3)
    
    pyplot.figure(figsize = (12, 10))
    pyplot.gca().set_ylim([1.0e1, 1.0e6])
    loglog(n_H, P)
    xlabel('log($n_H$)')
    ylabel('log(P)')
    pyplot.savefig(figure_name)
    print "\nPlot of pressure vs. density in thermal equilibrium was saved to: ", figure_name
    pyplot.close()

def log_equilibrium_temperature(n_H, cooling_function, heating_function, logT_min = 0.0, logT_max = 5.0):
    function = lambda logT: n_H * cooling_function(None, logT = logT) - heating_function()
    return bisect(function, logT_min, logT_max, 1.0e-8)

def bisect(function, a, b, tol):
    """
    Searches for a root between a and b using bisection, with an accuracy within tol.
    """
    fa = function(a)
    fb = function(b)
    if fa * fb > zero:
        raise Exception("Limits {0} and {1} may not enclose a root.".format(a, b))
    
    while (abs(a-b) > tol):
        c = 0.5*(a+b)
        fc = function(c)
        if fb * fc <= zero:
            a = c
            fa = fc
        else:
            b = c
            fb = fc
    
    if abs(fa) < abs(fb):
        return a
    else:
        return b


def mu(X = None, Y = 0.25, Z = 0.02, x_ion = 0.1):
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

global_mu = mu()

def u_from_T(T):
    return 3.0/2.0 * constants.kB * T / global_mu

def T_from_u(u):
    return 2.0/3.0 * u * global_mu / constants.kB

def evolve_and_plot_internal_energy(u_0, dt, n_H, T_0 = None, du_dt_adiabatic = zero, figure_name = None):
    if not T_0 is None:
        u_0 = u_from_T(T_0)
    function = lambda u: du_dt_adiabatic * u/u_0 + (gerritsen_heating_function() - n_H * my_cooling_function(T_from_u(u))) / global_mu
    time, u = integrate_ode_for_plot(function, u_0, dt)
    
    if figure_name and HAS_MATPLOTLIB:
        pyplot.figure(figsize = (12, 10))
        pyplot.subplot(2,1,1)
        plot(time, u.as_quantity_in(units.erg / units.g))
        xlabel('time')
        ylabel('internal energy')
        
        pyplot.subplot(2,1,2)
        plot(time, T_from_u(u))
        xlabel('time')
        ylabel('T')
        pyplot.savefig(figure_name)
        print "\nPlot of integrated ODE was saved to: ", figure_name
        pyplot.close()
    else:
        print u[-1]

def integrate_ode_for_plot(function, x, t_end, eps = 0.001):
    """
    Integrates the given ordinary differential equation of the form:
    dx/dt = function(x)
    for a time 't_end', using the initial value 'x'.
    The routine takes small steps, such that (abs(dx) <= eps * x)
    """
    t = 0 | units.s
    times = [0.0] | t_end.unit
    values = x.as_vector_with_length(1)
    while t < t_end:
        fx = function(x)
        step = min( (t_end-t), (eps*x)/abs(fx) )
        t += step
        x += fx * step
        times.append(t)
        values.append(x)
    return times, values



def evolve_internal_energy(u_0, dt, n_H, du_dt_adiabatic = zero):
    function = lambda u: (du_dt_adiabatic * u/u_0 + (gerritsen_heating_function() 
                        - n_H * my_cooling_function(T_from_u(u))) / global_mu)
    return integrate_ode(function, u_0, dt)

def integrate_ode(function, x, t_end, eps = 0.001):
    """
    Integrates the given ordinary differential equation of the form:
    dx/dt = function(x)
    for a time 't_end', using the initial value 'x'.
    The routine takes small steps, such that (abs(dx) <= eps * x)
    """
    t = 0 | units.s
    while t < t_end:
        fx = function(x)
        step = min( (t_end-t), ((eps*x)/abs(fx)).amin() )
        t += step
        x += fx * step
    return x



if __name__ == '__main__':
    set_printing_strategy("cgs")
    plot_cooling_function(gerritsen_cooling_function, "gerritsen_cooling_function.png")
    plot_cooling_function(my_cooling_function, "gerritsen_cooling_function_fudged.png")
    
    plot_cooling_vs_heating_n(0.1 | units.cm**-3, gerritsen_cooling_function, gerritsen_heating_function, "gerritsen_cooling_vs_heating_n-1.png")
    print "Equilibrium temperature:",
    print 10.0**log_equilibrium_temperature(0.1 | units.cm**-3, gerritsen_cooling_function, gerritsen_heating_function) | units.K
    
    plot_cooling_vs_heating_P(1.0e3 | units.K * units.cm**-3, my_cooling_function, gerritsen_heating_function, "gerritsen_cooling_vs_heating_P3.png")
    plot_cooling_vs_heating_P(1.0e2 | units.K * units.cm**-3, my_cooling_function, gerritsen_heating_function, "gerritsen_cooling_vs_heating_P2.png")
    plot_cooling_vs_heating_P(1.0e4 | units.K * units.cm**-3, my_cooling_function, gerritsen_heating_function, "gerritsen_cooling_vs_heating_P4.png")
    
    plot_pressure_density_in_thermal_equilibrium(gerritsen_cooling_function, gerritsen_heating_function, "pressure_density_gerritsen_cooling.png")
    plot_pressure_density_in_thermal_equilibrium(my_cooling_function, gerritsen_heating_function, "pressure_density_my_cooling.png")
    
    for i in range(-2, 3):
        dens = 10.0**i | units.cm**-3
        evolve_and_plot_internal_energy(None, 1e6 | units.yr, dens, T_0 = 1.0 | units.K, 
            du_dt_adiabatic = zero, figure_name = "internal_energy_evolution_{0:=03}.png".format(i))
        print "Equilibrium temperature:",
        print 10.0**log_equilibrium_temperature(dens, my_cooling_function, gerritsen_heating_function) | units.K
    
    for i in range(-2, 3):
        dens = 10.0**i | units.cm**-3
        evolve_and_plot_internal_energy(None, 1.0 | units.Myr, dens, T_0 = 1.0 | units.K, 
            du_dt_adiabatic = u_from_T(10 | units.K)**1.5 / (1.0 | units.parsec), 
            figure_name = "internal_energy_evolution_dudtad_{0:=03}.png".format(i))
        print "Equilibrium temperature:",
        print 10.0**log_equilibrium_temperature(dens, my_cooling_function, gerritsen_heating_function) | units.K
    
    evolve_and_plot_internal_energy(None, 10.0 | units.Myr, 1.0 | units.cm**-3, T_0 = 1.0 | units.K, 
        du_dt_adiabatic = u_from_T(10 | units.K)**1.5 / (1.0 | units.parsec), 
        figure_name = "internal_energy_evolution_dudtad_000_long.png")
    
    u_0 = u_from_T([1000.0]*5 | units.K)
    n_H = (units.cm**-3).new_quantity([10.0**i for i in range(-2, 3)])
    print T_from_u(evolve_internal_energy(u_0, 1e6 | units.yr, n_H))
    
