#!/usr/bin/env python

"""
This module contains the temperature and pressure profile
of the 1976 US standard atmosphere model up until 47,000m.

Reference
---------
http://en.wikipedia.org/wiki/US_Standard_Atmosphere

Parameters
----------
P0 : pressure at sea level
T0 : temperature at sea level
B1 : Temperature Lapse Rate between 0m and H1
     -6.5 K/km
B2 : Temperature Lapse Rate between H2 and H3
      1 K/km
B3 : Temperature Lapse Rate beyond H3
      2.8 K/km
H1 : altitude of first ladder (~11,000m)
H2 : altitude of second ladder (~20,000m)
H3 : altitude of third ladder (~32,000m)
"""

# Parameters of the model (pressures in Pa
# temperatures in K)
# P0 and T0 are taken to be the ones corresponding
# to the actual atmospheric conditions on the day
# of the balloon flight
P0, T0 = 103200, 277.444
H1, H2, H3 = 11000, 20000, 32000
B1, B2, B3 = -0.0065, 0.001, 0.0028


def temperature(z):
    """
    Parameters
    ----------
    z: float
        altitude in meters (m)

    Returns
    -------
    temp: float
        temperature in Kelvin (K)
    """
    if abs(z) <= H1:
        return T0 + B1 * z
    elif H1 < abs(z) <= H2:
        return T0 + B1 * H1
    elif H2 < abs(z) <= H3:
        return T0 + B1 * H1 + B2 * (abs(z) - H2)
    else:
        T2 = T0 + B1 * H1 + B2 * (H3 - H2)
        return T2 + B3 * (abs(z) - H3)


def pressure(z):
    """
    The pressure is obtained by integrating the hydrostatic
    equation indicating local mechanical equilibrium of the
    atmosphere, along with the use of the equation of state
    for Perfect Gases.

    dP
    -- = - rho * g,            P = rho * R * T / Mair
    dz

    P is the pressure
    z the altitude
    rho is the volumic mass of the air
    g is the acceleration of gravity (~9.807 m/s**2)
    R is the Perfect Gas constant (~8.31446 J/(mol.K))
    T is the temperature
    Mair is the molar mass of the air (~28.97 g/mol)
    """

    from math import exp

    # convenient constant a = (Mair * g) / R
    a = 0.034170444021620165

    if abs(z) <= H1:
        return P0 * (1. + B1 * abs(z) / T0) ** (-a / B1)
    elif H1 < abs(z) <= H2:
        P1 = P0 * (1. + B1 * H1 / T0) ** (-a / B1)
        T1 = temperature(H1)
        return P1 * exp(-a * (abs(z) - H1) / T1)
    elif H2 < abs(z) <= H3:
        P1 = P0 * (1. + B1 * H1 / T0) ** (-a / B1)
        T1 = temperature(H1)
        P2 = P1 * exp(-a * (H2 - H1) / T1)
        return P2 * (1. + B2 * (abs(z) - H2) / T1) ** (-a / B2)
    else:
        P1 = P0 * (1. + B1 * H1 / T0) ** (-a / B1)
        T1 = temperature(H1)
        P2 = P1 * exp(-a * (H2 - H1) / T1)
        T2 = temperature(H3)
        P3 = P2 * (1. + B2 * (H3 - H2) / T1) ** (-a / B2)
        return P3 * (1. + B3 * (abs(z) - H3) / T2) ** (-a / B3)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    z_list = np.linspace(0.,40000.,40001)
    temp = np.array([temperature(z) for z in z_list])
    pres = np.array([pressure(z) for z in z_list])
    fig = plt.figure()
    ax = fig.add_subplot(2, 1, 1)
    ax.set_ylabel(r'altitude z (m)', fontsize=14)
    ax.set_xlabel(r'temperature (K)', fontsize=12)
    ax.plot(temp, z_list, '-r', label='temperature profile')

    ax = fig.add_subplot(2, 1, 2)
    ax.set_ylabel(r'altitude z (m)', fontsize=14)
    ax.set_xlabel(r'pressure (Pa)', fontsize=12)
    plt.plot(pres, z_list, '-b', label='pressure profile')

    plt.show()
