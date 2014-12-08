#!/usr/bin/env python

"""
Uses the 'curvature flow' - a 1d example of a Ricci flow - in order to
evolve the shape of a 1d closed curve in the cartesian plane, so as to
collapse it to a point.

The main feature of the 'curvature flow' is the finite time to collapse
along with nice properties like the 'rounding up' of the curve as time
evolves.

The master equation is exactly a heat equation. For each point on the
curve at time t, with coordinates `x(s, t)` and `y(s, t)` where `s` is
the arclength coordinate along the curve, this equation reads:
    D(z(s, t), t) = D(z(s, t), s, 2)
where `z=(x,y)`.

We chose a discretization scheme in which the time coordinate `t` is
treated with a simple Euler first order scheme, and the space coordinate
`s`, is treated with a centered second-order scheme.
"""

from math import sqrt
import numpy as np
from scipy.interpolate import interp1d
from scipy.sparse import spdiags


def curve_flow(x0, y0, dt=0.01, nstep=100):

    # -- conservative choice for the length step
    #    and number of time step given the input
    #    time step. With this length step we make
    #    sure that the algorithm is stable
    ds = sqrt(2. * dt)

    x0 = np.array(x0).astype('f')
    y0 = np.array(y0).astype('f')

    assert x0.ndim == 1
    assert y0.ndim == 1
    assert x0.size == y0.size

    assert 0. < dt
    assert 0. < ds

    x_coords = []
    y_coords = []

    xold, yold = get_equispaced_pts(x0, y0, ds)

    for n in xrange(nstep):
        print 'step %6i' % n
        if xold.size <= 10:
            break
        else:
            x_coords += [xold]
            y_coords += [yold]
            xnew, ynew = flow(xold, yold, dt, ds)
            xold, yold = get_equispaced_pts(xnew, ynew, ds)

    return x_coords, y_coords


def flow(xold, yold, dt, ds):

    itemtype = xold.dtype
    n = xold.size

    # -- non-empty diagonals of the matrix
    near = (1. * dt / ds**2) * np.ones(n, dtype=itemtype)
    diag = (1. - 2. * dt / ds**2) * np.ones(n, dtype=itemtype)
    data = np.vstack((diag, near, near, near, near))

    # -- diagonal indices of the non-empty diagonals
    indices = np.array([0, -1, 1, -(n - 1), n - 1])

    A = spdiags(data, indices, n, n)

    xnew = A.dot(xold)
    ynew = A.dot(yold)

    return xnew, ynew


def get_equispaced_pts(xold, yold, ds):

    itemtype = xold.dtype

    # -- old grid of arclength values
    sold = get_arclength_grid(xold, yold)
    xxold = np.append(xold, xold[0])
    yyold = np.append(yold, yold[0])

    fx = interp1d(sold, xxold, kind='cubic')
    fy = interp1d(sold, yyold, kind='cubic')

    snew = np.arange(sold.min(), sold.max(), ds, dtype=itemtype)

    xnew = fx(snew)
    ynew = fy(snew)

    return xnew, ynew


def get_arclength_grid(x, y):

    assert x.ndim == 1
    assert y.ndim == 1
    assert x.size == y.size

    s_out = np.zeros(x.size + 1, dtype=x.dtype)

    xr = np.roll(x, -1)
    yr = np.roll(y, -1)

    s_int = np.sqrt((xr - x)**2 + (yr - y)**2)

    s_out[1:] = np.cumsum(s_int)

    return s_out


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    from math import pi

    # -- test on a circle
    R = 0.5
    theta = np.arange(0., 2. * pi, 0.03)
    x0 = R * np.cos(theta)
    y0 = R * np.sin(theta)

    x_coords, y_coords = curve_flow(x0, y0, dt=0.001)
    center = len(x_coords) / 2

    plt.plot(x_coords[0], y_coords[0], 'b',
             x_coords[center], y_coords[center], 'g',
             x_coords[-1], y_coords[-1], 'r')
    plt.show()
