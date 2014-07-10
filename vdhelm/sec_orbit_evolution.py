"""
    Formulas to calculate the secular orbital evolution.
    Mass loss formulas based on Sepinsky, Willems Kalogera and Rasio 2007 and 2009
    Tidal formulas based on Hurley et al. 2002 after Hut 1981, Zahn 1977 and Rasio 1996
"""
import numpy

from amuse.support.console import set_printing_strategy
from amuse.units import units, constants
from amuse.ext.roche_radius import Roche_Orbit

class Evolving_Orbit(Roche_Orbit):
    def __init__(self,
                mass_loss=0.|units.MSun/units.yr,
                cos_accretion_angle_P=1.,
                radius_1=None,
                radius_2=1|units.RSun,
                accreted_fraction=1.,
                lost_angular_momentum=None,
                envelope_mass_fraction=0.8,
                envelope_radius_fraction=0.95,
                radius_of_gyration_squared=0.1,
                is_convective=True,
                luminosity=1|units.LSun,
                **kwargs):
        super(Evolving_Orbit, self).__init__(**kwargs)
        self.mass_loss = mass_loss
        self.cos_accretion_angle_P = cos_accretion_angle_P
        self._radius_1 = radius_1
        self.radius_2 = radius_2
        self.accreted_fraction = accreted_fraction
        self._lost_angular_momentum = lost_angular_momentum
        self.envelope_mass_fraction = envelope_mass_fraction
        self.envelope_radius_fraction = envelope_radius_fraction
        self.radius_of_gyration_squared = radius_of_gyration_squared # See Hut (1981) eq. A26
        self.is_convective = is_convective
        self.luminosity = luminosity

    @property
    def lost_angular_momentum(self):
        """
            The mu parameter is the angular momentum lost through mass loss from the system.
            If lost_angular_momentum is set to None, this is M1/M2 following Kolb et al. 2001.
        """
        if self._lost_angular_momentum is None:
            return self.mass_ratio
        else:
            return self._lost_angular_momentum

    @lost_angular_momentum.setter
    def lost_angular_momentum(self, value):
        self._lost_angular_momentum = value

    @property
    def radius_1(self):
        """
            Radius of the donor star, return the Roche radius if it is set to None
        """
        if self._radius_1 is None:
            return self.sepinsky_roche_radius()
        else:
            return self._radius_1

    @radius_1.setter
    def radius_1(self, value):
        self._radius_1 = value

    def envelope_mass(self):
        return self.mass_1 * self.envelope_mass_fraction

    def envelope_radius(self):
        return self.radius_1 * self.envelope_radius_fraction

    def orbital_angular_velocity(self):
        """ From Sepinsky 2007 """
        return 2 * numpy.pi / self.period * (1 + self.eccentricity * numpy.cos(self.true_anomaly))**2 / (1-self.eccentricity**2)**(3./2.)

    def mean_orbital_angular_velocity(self):
        """ From Hut 1981 """
        return (constants.G * (self.mass_1 + self.mass_2))**(1./2.) * self.semimajor_axis**(-3./2.)

    def peri_orbital_angular_velocity(self):
        """ From Hut 1981 """
        return (1 + self.eccentricity)**2 * (1-self.eccentricity**2)**(-3./2.) * self.mean_orbital_angular_velocity()

    def angular_velocity_over_mean_orbital(self):
        """ Note that Sepinsky and Hut use different angular velocity ratio definitions """
        return self.angular_velocity_ratio * self.peri_orbital_angular_velocity() / self.mean_orbital_angular_velocity()

    def mass_loss_semimajor_axis_rate(self):
        """ Sepinsky Equation 39 / Equation 18 """
        num_0 = self.semimajor_axis * self.mass_loss
        den_0 = numpy.pi * self.mass_1 * numpy.sqrt(1 - self.eccentricity**2)

        part_1 = self.eccentricity * self.L1_radius() / self.semimajor_axis

        part_2 = self.accreted_fraction * self.mass_ratio * self.eccentricity * self.radius_2 * self.cos_accretion_angle_P / self.semimajor_axis

        part_3 = (self.accreted_fraction * self.mass_ratio - 1) * (1 - self.eccentricity**2)

        part_4 = (1. - self.accreted_fraction) * (self.lost_angular_momentum + 0.5) * (1 - self.eccentricity**2) * self.mass_ratio / (1 + self.mass_ratio)

        return num_0 / den_0 * ( part_1 + part_2 + part_3 + part_4)

    def mass_loss_eccentricity_rate(self):
        """ Sepinsky Equation 40 """
        num_0 = numpy.sqrt(1 - self.eccentricity**2) * self.mass_loss
        den_0 = 2. * numpy.pi * self.mass_1

        part_1 = self.accreted_fraction * self.mass_ratio * self.radius_2 * self.cos_accretion_angle_P / self.semimajor_axis

        part_2 = self.L1_radius() / self.semimajor_axis

        part_3 = 2. * (self.accreted_fraction * self.mass_ratio - 1) * (1 - self.eccentricity)

        part_4 = 2. * (1. - self.accreted_fraction) * (self.lost_angular_momentum + 0.5) * (1 - self.eccentricity) * self.mass_ratio / (1 + self.mass_ratio)

        return num_0 / den_0 * ( part_1 + part_2 + part_3 + part_4)

    def tidal_semimajor_axis_rate(self):
        """ Hut eq 9 """
        part_1 = -6 * self.k_over_T() * self.q_2() * (1 + self.q_2())
        part_2 = (self.radius_1 / self.semimajor_axis) ** 8 * self.semimajor_axis / ((1 - self.eccentricity**2)**7.5)
        part_3 = self.f_1() - (1 - self.eccentricity**2)**1.5 * self.f_2() * self.angular_velocity_over_mean_orbital()

        return part_1 * part_2 * part_3

    def tidal_eccentricity_rate(self):
        """ Hurley equation 25 after Hut eq 10 """
        part_1 = -27 * self.k_over_T() * self.q_2() * (1 + self.q_2())
        part_2 = (self.radius_1 / self.semimajor_axis) ** 8 * self.eccentricity / ((1 - self.eccentricity**2)**6.5)
        part_3 = self.f_3() - 11./18. * (1 - self.eccentricity**2)**1.5 * self.f_4() * self.angular_velocity_over_mean_orbital()

        return part_1 * part_2 * part_3

    def tidal_rotation_rate(self):
        """ Hurley equation 26 after Hut eq 11"""
        part_1 = 3 * self.k_over_T() * self.q_2()**2 / self.radius_of_gyration_squared
        part_2 = (self.radius_1 / self.semimajor_axis) ** 6 * self.mean_orbital_angular_velocity() / ((1 - self.eccentricity**2)**6)
        part_3 = self.f_2() - (1 - self.eccentricity**2)**1.5 * self.f_5() * self.angular_velocity_over_mean_orbital()

        return part_1 * part_2 * part_3

    def k_over_T(self):
        if self.is_convective:
            return self.convective_k_over_T()
        else:
            return self.radiative_k_over_T()

    def convective_k_over_T(self):
        """ Hurley equation 30 after Rasio 1996 """
        f_conv = 1 # This numerical factor depends on the tides but it is usually 1.
        eddy_turnover_timescale = 0.4311 * ( self.envelope_mass() * self.envelope_radius() * ( self.radius_1 - 0.5 * self.envelope_radius()) / (3 * self.luminosity) )**(1./3.)

        return 2. / 21. * f_conv / eddy_turnover_timescale * self.envelope_mass() / self.mass_1

    def radiative_k_over_T(self):
        """ Hurley equation 42 after Zahn 1975 """
        # The units are not explicitly removed in Hurley, but this seems to be the only way to match the numbers from Zahn 1975 this is based on.
        M = self.mass_1.value_in(units.MSun)
        R = self.radius_1.value_in(units.RSun)
        a = self.semimajor_axis.value_in(units.RSun)
        tidal_coefficient = 1.592e-9 * M**2.84
        return (1.9782e4 * M * R**2 / a**5 * (1 + self.q_2())**(5./6.) * tidal_coefficient) | units.yr**-1

    def q_2(self):
        return self.mass_2 / self.mass_1

    def alpha(self):
        """ Based on Hut 1981 equation 22 """
        #TODO: how do you calculate this a_0?
        equilibrium_semimajor_axis = self.semimajor_axis
        return self.q_2()/(1.+self.q_2()) * 1./self.radius_of_gyration_squared * (equilibrium_semimajor_axis/self.radius_1)**2

    def f_1(self):
        e = self.eccentricity
        return 1 + 31/2*e**2 + 255/8*e**4 + 185/16*e**6 + 25/64*e**8

    def f_2(self):
        e = self.eccentricity
        return 1 + 15/2*e**2 + 45/8*e**4 + 5/16*e**6

    def f_3(self):
        e = self.eccentricity
        return 1 + 15/4*e**2 + 15/8*e**4 + 5/64*e**6

    def f_4(self):
        e = self.eccentricity
        return 1 + 3/2*e**2 + 1/8*e**4

    def f_5(self):
        e = self.eccentricity
        return 1 + 3*e**2 + 3/8*e**4

    def L1_radius(self):
        """ Equation A15 Sepinsky 2007b """
        log_q = numpy.log(self.mass_ratio)
        f = self.angular_velocity_ratio
        e = self.eccentricity
        X_L1 = 0.529 + 0.231 * log_q - f**2 * ( 0.031 + 0.025 * e) * (1 + 0.4 * log_q)
        return self.separation() * X_L1

if __name__ == '__main__':
    set_printing_strategy("custom", preferred_units = [units.MSun, units.RSun, units.day, units.LSun, units.AU *units.day**-1, units.day**-1 * units.day**-1])

    simple_test()
