#!/usr/bin/env python

# Authors: Nicolas Poilvert
# License: BSD

"""
The purpose of this program is to be able to compute "lattice sums" which first
appeared in an old paper by Lord Rayleigh. This gentleman had the brilliant
insight to come up with what became to be known as the "Rayleigh Identity",
which introduced the notion of "lattice sums".

Reference papers
================

Lord Rayleigh, Phil. Mag., 34, 481-502 (1892)

H. Cheng and S. Torquato, Mathematical, Physical and Engineering Sciences,
                          453, 145-161 (1997)
"""

import numpy as np
from scipy.misc import factorial2
from scipy.special import gammainc

# -- Defining a special simplified kind of spherical harmonics (defined in
# Greengard's thesis from MIT)
def Ylm(l, m, theta, phi):
    """
    The definition we use in our code is the one of Greengard. Namely, we have
    the following expression for Ylm(t,p)::
                                               |m|
        Ylm(t, p) = sqrt((l-|m|)! / (l+|m|)!).P  (cos(t)).exp(i.m.p)
                                               l
    The spherical harmonics defined in scipy is (scipy calls "theta" what we
    consider to be phi and vice-versa)::
         S
        Ylm(p, t) = sqrt(2l+1 / 4.pi).Ylm(t, p) if m >= 0
         S
        Ylm(p, t) = (-1)**m.sqrt(2l+1 / 4.pi).Ylm(t, p) if m < 0
    """

    from scipy.special import sph_harm

    assert l >= 0
    assert abs(m) <= l

    prefactor = np.sqrt(4.*np.pi / (2.*l+1.))

    if m >= 0:
        return prefactor * sph_harm(m, l, phi, theta)
    else:
        return prefactor * np.conjugate(sph_harm(abs(m), l, phi, theta))

# -- some necessary utilities to deal with the direct and reciprocal spaces of
# cubic crystals (Simple Cubic, Body-Centered Cubic and Face-Centered Cubic)
def direct_space_elementary_vectors(a=1., crystal_type='sc'):
    """
    This routine will generate the elementary basis vectors (a1, a2, a3) for a
    given cubic crystal Bravais lattice type ('simple cubic', 'body-centered
    cubic', 'face-centered cubic'
    """

    assert a > 0.
    assert crystal_type in ['sc', 'bcc', 'fcc']

    if crystal_type == 'sc':
        a1 = a * np.array([1., 0., 0.])
        a2 = a * np.array([0., 1., 0.])
        a3 = a * np.array([0., 0., 1.])
    elif crystal_type == 'bcc':
        a1 = a/2. * np.array([1., 1., -1.])
        a2 = a/2. * np.array([-1., 1., 1.])
        a3 = a/2. * np.array([1., -1., 1.])
    elif crystal_type == 'fcc':
        a1 = a/2. * np.array([1., 1., 0.])
        a2 = a/2. * np.array([0., 1., 1.])
        a3 = a/2. * np.array([1., 0., 1.])
    else:
        raise ValueError('unsupported crystal type')

    return a1, a2, a3

def Cartesian2Spherical(cartesian_coordinates):

    assert cartesian_coordinates.ndim == 2
    N, dim = cartesian_coordinates.shape
    assert dim == 3

    x, y, z = cartesian_coordinates.T

    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    elevation = np.arctan2(hxy, z)
    longitude = np.arctan2(y, x)

    spherical_coordinates = np.vstack([r, elevation, longitude]).T

    return spherical_coordinates

def Spherical2Cartesian(spherical_coordinates):

    assert spherical_coordinates.ndim == 2
    N, dim = spherical_coordinates.shape
    assert dim == 3

    r, theta, phi = spherical_coordinates.T

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    cartesian_coordinates = np.vstack([x, y, z]).T

    return cartesian_coordinates

def lattice_point_spherical_coordinates(a=1., crystal_type='sc', N=10):
    """
    Given a crystal type and the lattice spacing 'a', and given a size for the
    lattice of points about the origin (same size in each directions x, y and
    z), this function returns the spherical coordinates of all the lattice
    points (except the origin).

    The origin is assumed to be (0., 0., 0.).
    """

    a1, a2, a3 = direct_space_elementary_vectors(a=a, crystal_type=crystal_type)

    coordinates = []
    for i in xrange(-N, N+1):
        for j in xrange(-N, N+1):
            for k in xrange(-N, N+1):
                if (i, j, k) != (0, 0, 0):
                    coordinates += [i*a1 + j*a2 + k*a3]
    coordinates = np.array(coordinates)
    coordinates = Cartesian2Spherical(coordinates)
    assert coordinates.shape == ((2*N+1)**3-1, 3)

    return coordinates

# -- Core routine to compupte lattice sums from direct summations
def lattice_sum(l, m, crystal_type='sc', N=10):
    """
    Computes the necessary lattice sums for cubic crystals::

        Sum[ Ylm(theta_j, phi_j) / r_j ** (l+1), j=-1,1,-2,2,...]

    except if l = 2 and m = 0 for which we must use the physically relevant
    conditionally converging values::

        2.pi/3 for SC
        4.pi/3 for BCC
        8.pi/3 for FCC
    """
    assert l >= 2
    assert crystal_type in ['sc', 'bcc', 'fcc']

    if l == 2 and m == 0:
        if crystal_type == 'sc':
            return 2.*np.pi/3.
        elif crystal_type == 'bcc':
            return 4.*np.pi/3.
        elif crystal_type == 'fcc':
            return 8.*np.pi/3.
        else:
            raise ValueError('unsupported crystal type')
    else:
        spherical_coordinates = lattice_point_spherical_coordinates(
                    a=1., crystal_type=crystal_type, N=N
                )
        r, theta, phi = spherical_coordinates.T
        radial_values = 1. / r ** (l+1)
        angular_values = Ylm(l, m, theta, phi)
        lattice_sum_value = np.sum(radial_values * angular_values)
        return lattice_sum_value

# -- A much smarter way to compute the sum using Ewald's summations
def Ewald_lattice_sum(l, m, crystal_type='sc', N=5, alpha=np.pi):
    """
    Accelerated lattice sums using the Ewald method from S. Adler, Physica 27,
    pp 1193-1201 (1961), equation (25) on page 1200.
    """
    assert l >= 2
    assert crystal_type in ['sc', 'bcc', 'fcc']

    if l == 2 and m == 0:
        if crystal_type == 'sc':
            return 2.*np.pi/3.
        elif crystal_type == 'bcc':
            return 4.*np.pi/3.
        elif crystal_type == 'fcc':
            return 8.*np.pi/3.
        else:
            raise ValueError('unsupported crystal type')
    else:
        # -- Spherical coordinates of the direct-space lattice points
        spherical_coordinatesj = lattice_point_spherical_coordinates(
                    a=1., crystal_type=crystal_type, N=N
                )
        rj, thetaj, phij = spherical_coordinatesj.T
        angular_valuesj = Ylm(l, m, thetaj, phij)

        # -- Volume of the unit cell (used later on)
        if crystal_type == 'sc':
            V = 1.
        elif crystal_type == 'bcc':
            V = 1. / 2.
        elif crystal_type == 'fcc':
            V = 1. / 4.

        # -- Spherical coordinates of the reciprocal-space lattice points
        if crystal_type == 'bcc':
            spherical_coordinatesn = lattice_point_spherical_coordinates(
                        a=4.*np.pi, crystal_type='fcc', N=N
                    )
        elif crystal_type == 'fcc':
            spherical_coordinatesn = lattice_point_spherical_coordinates(
                        a=4.*np.pi, crystal_type='bcc', N=N
                    )
        elif crystal_type == 'sc':
            spherical_coordinatesn = lattice_point_spherical_coordinates(
                        a=2.*np.pi, crystal_type='sc', N=N
                    )
        else:
            raise ValueError('unsupported crystal type')
        gn, thetan, phin = spherical_coordinatesn.T
        angular_valuesn = Ylm(l, m, thetan, phin)

        # -- now computing the lattice sum in two terms. The first is a sum over
        # the direct-space lattice while the second is a sum over the
        # reciprocal-space lattice
        lattice_sumj = np.sum(
                    angular_valuesj * (1. / rj) ** (l+1) * \
                    (1. - gammainc(l+0.5, alpha*rj**2))
                )
        lattice_sumn = (1j)**l * 4.*np.pi / V * np.sum(
                    angular_valuesn * gn**(l-2) * (1./factorial2(2*l-1,
                        exact=False)) * np.exp(-gn**2/(4.*alpha))
                )
        lattice_sum_value = lattice_sumj + lattice_sumn

        return lattice_sum_value
