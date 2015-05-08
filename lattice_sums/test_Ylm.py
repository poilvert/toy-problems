#!/usr/bin/env python

# Authors: Nicolas Poilvert
# License: BSD

import numpy as np
from cubic_3D_lattices import Ylm

# ------------------------------------------------------------------------------
# -- utility to generate appropriate theta, phi arrays to test orthonormality of
# the spherical harmonics
# ------------------------------------------------------------------------------
def theta_phi_arrays(N):

    # -- 1d arrays for theta and phi angles
    theta_1d_array = np.linspace(0., np.pi, num=(N+1))
    phi_1d_array = np.linspace(0., 2.*np.pi, num=(2*N+1))
    # -- angle increments
    d_theta = np.pi/N
    d_phi = 2.*np.pi/(2.*N)
    # -- 2d arrays for theta and phi angles (theta, phi) tables
    phi_2d_array, theta_2d_array = np.meshgrid(phi_1d_array, theta_1d_array)
    # -- sinus of theta angle for the 2d case (used later on for computing the
    # solid angle element)
    sin_theta = np.sin(theta_2d_array)
    # -- solid angle infinitesimal elements (for the 2d arrays)
    solid_angle_elements = d_theta * d_phi * sin_theta

    return theta_2d_array, phi_2d_array, solid_angle_elements


def inner_product(l1, m1, l2, m2, N=20):

    # -- prepare the necessary (theta, phi) arrays
    theta2d, phi2d, dOmega = theta_phi_arrays(N)

    # -- compute spherical harmonics onto those arrays
    Ylm1_values = Ylm(l1, m1, theta2d, phi2d)
    Ylm2_values = Ylm(l2, m2, theta2d, phi2d)

    # -- integration
    integral = np.sum(np.conjugate(Ylm1_values) * Ylm2_values * dOmega)

    return integral

# --------
# -- Tests
# --------
def test_same_spherical_harmonics_m_negative():

    N = 200
    l, m = 2, -1
    norm_squared = inner_product(l, m, l, m, N=N)
    expected_norm_squared = 4. * np.pi / (2.*l+1.)

    assert np.abs(norm_squared - expected_norm_squared) < 1e-2

def test_same_spherical_harmonics_m_positive():

    N = 100
    l, m = 3, 2
    norm_squared = inner_product(l, m, l, m, N=N)
    expected_norm_squared = 4. * np.pi / (2.*l+1.)

    assert np.abs(norm_squared - expected_norm_squared) < 1e-2

def test_different_spherical_harmonics_same_l():

    N = 100
    l1, m1 = 3, 2
    l2, m2 = 3, -1
    inner_prod = inner_product(l1, m1, l2, m2, N=N)
    expected_inner_prod = 0.

    assert np.abs(inner_prod - expected_inner_prod) < 1e-2

def test_different_spherical_harmonics_same_m():

    N = 100
    l1, m1 = 1, 1
    l2, m2 = 3, 1
    inner_prod = inner_product(l1, m1, l2, m2, N=N)
    expected_inner_prod = 0.

    assert np.abs(inner_prod - expected_inner_prod) < 1e-2
