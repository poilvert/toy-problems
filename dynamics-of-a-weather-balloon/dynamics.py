#!/usr/bin/env python

"""
Module to solve for the equation of motion of the balloon
as it ascends.

Assumptions:
    1. motion only along z (lateral forces completely
    decoupled from the vertical motion)
    2. No mass loss (i.e. no loss of Helium through
    the balloon membrane)
    3. Spherical approximation for the shape of the
    balloon at every altitude.
    4. External forces are "weight" and "buoyancy" from
    the air and air drag (quadratic in speed)

Models:
    1. An 'adiabatic' model in which the motion of the balloon
    is so fast (empirically it was an estimate ~ 0.5h between
    ascent and descent so quite fast!) that the temperature
    of the Helium inside the balloon does not change. The
    constant temperature is denoted T0.
    2. A 'quasi-static' model in which the temperature of the
    Helium gas inside the balloon always equilibrates with the
    temperature of the outside air at every altitude.

Parameters:
    - M0 mass of balloon (balloon + equipment)
    - M mass of Helium in the balloon
    - Cd drag coefficiant (for a sphere at not too
      high and not too low Reynolds number ~ 0.5)
    - Mair molecular weight of air (kg/mol)
    - MHe molecular weight of Helium (kg/mol)
    - R perfect gas constant
    - g acceleration of gravity near the surface of the
      earth
"""

from math import pi
import numpy as np
from scipy.integrate import odeint

# getting the experimental data
from parse_raw_data import parse_data

ascent_alt, descent_alt, ascent_time, descent_time = parse_data()

# importing Temperature and Pressure functions
from atmosphere import temperature as Tair
from atmosphere import pressure as Pair


# parameters of the balloon (SI units)
M = 0.75       # mass of Helium in baloon (kg)
g = 9.807      # acceleration of gravity
Mair = 0.02897 # molecular weight of air
MHe = 0.004    # molecular weight of Helium
R = 8.31446    # perfect gas constant

# (T0, M0, Cd), parameters of the model to fit
# T0 Helium temperature in Kelvins
# M0 empty balloon mass in kg
# Cd drag coefficiant (no dimension)
T0_range = np.arange(277, 278, 1)
M0_range = np.arange(0.1, 1.1, 0.1)
Cd_range = np.arange(0.4, 1.25, 0.05)

# time range for the ODE (in seconds)
tmin = ascent_time.min()
tmax = ascent_time.max()
time = np.arange(tmin, tmax+1, 0.1)

# solver
def solver(t_range, T0, M0, Cd, model='adiabatic'):
    """
    this solver takes a set of parameters (T0, M0, Cd) and
    returns the integrated solution z(t) from the equation
    of motion.

    Parameters
    ----------
    t_range: 1D array of floats
        time interval over which to compute the solution
    T0: 1D array of integers
        grid of temperature onto which the fit takes place
    M0: 1D array of floats
        grid of load mass onto which to perform the fit
    Cd: 1D array of floats
        grid of drag coefficiant values onto which to perform
        the fit
    model: string
        which type of model to use for the equation of motion

    Returns
    -------
    X: list of 1D arrays
        [z(t), z'(t)], respectively the altitude and velocity
        functions with time
    """
    # initial conditions (no altitude, no velocity)
    X0 = [0., 0.]

    # this function gives the gradient of the dynamical
    # system (one needs to transform the equation of motion which is
    # second order, into a system of first order ODEs) to use the
    # generic 'odeint' integrator from scipy.
    # if X(t) = [z(t), z'(t)], then
    # X'(t) = gradient(X(t)), where the 'gradient' function is obtained
    # from the equation of motion.
    if model == 'adiabatic':
        def gradient(X, t):
            alpha = M/(M0+M)*(Mair/MHe)*g*T0
            beta = 0.5*Cd/(M0+M)*(Mair/R)*pi*(3./(4.*pi)*(M/MHe)*R*T0)**(2./3.)
            return [X[1],
                    -g+alpha/Tair(X[0])-beta*(Pair(X[0])**(1./3.)/Tair(X[0]))*(X[1])**2]
    elif model == 'quasi-static':
        def gradient(X, t):
            alpha = M/(M0+M)*(Mair/MHe)*g
            beta = 0.5*Cd/(M0+M)*(Mair/R)*pi*(3./(4.*pi)*(M/MHe)*R)**(2./3.)
            return [X[1],
                    -g+alpha-beta*(Pair(X[0])/Tair(X[0]))**(1./3.)*(X[1])**2]

    # solution
    X = odeint(gradient, X0, t_range)

    return X

# time values for which to compute the solution corresponds to
# the original experimental values
time_eval = []
for t in time:
    if t in ascent_time:
        time_eval += [True]
    else:
        time_eval += [False]
time_eval = np.array(time_eval, dtype=bool)

assert time_eval.size == time.size
assert np.unique(ascent_time).size == ascent_time.size
assert time_eval.sum() == ascent_time.size

# computing the standard deviation between the theoretically predicted
# altitude values and the experimental ones
error = np.zeros((len(T0_range), len(M0_range), len(Cd_range)),
                 dtype=np.float64)
for i, T0 in enumerate(T0_range):
    for j, M0 in enumerate(M0_range):
        for k, Cd in enumerate(Cd_range):
            X = solver(time, T0, M0, Cd, model='quasi-static')
            alt, vel = X.T
            solution = alt[time_eval]
            assert solution.size == ascent_alt.size
            error[i, j, k] = (solution - ascent_alt).std()

solution_argmin = error.argmin()
optimal_indices = np.unravel_index(solution_argmin, error.shape)

Topt = T0_range[optimal_indices[0]]
Mopt = M0_range[optimal_indices[1]]
Copt = Cd_range[optimal_indices[2]]

print 'optimal solution :'
print ' - T0 = %5i K' % Topt
print ' - M0 = %5.2f kg' % Mopt
print ' - Cd = %5.2f' % Copt
print 'root mean square deviation from experiment :'
print ' - std = %10.2f' % error[optimal_indices[0], optimal_indices[1],
                                optimal_indices[2]]

import matplotlib.pyplot as plt
plt.plot(ascent_time, ascent_alt, 'r-', label='experiment')
X = solver(time, Topt, Mopt, Copt)
alt, vel = X.T
solution = alt[time_eval]
plt.plot(ascent_time, solution, 'b-', label='theory (best match)')
plt.xlabel('time (s)')
plt.ylabel('altitude (m)')
plt.grid(True)
plt.legend(loc=2)
plt.show()
