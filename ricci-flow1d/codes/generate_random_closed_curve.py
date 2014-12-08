#!/usr/bin/env python

from math import pi
import numpy as np


def random_closed_curve(theta_step=0.01, a=1., b=1.2, c= 0.5, d=0.6):
    """A simple program that generates a simple closed planar curve
    using polar coordinates. The program returns the x and y coordinates
    of points along the curve. The points are regularly spaced in *angle*
    which does not necessarily means that they are regularly spaced
    along the curve (i.e. in terms of arclength).
    """

    def radius_fct(theta_range, a, b, c, d):
        res = a * np.cos(theta_range) + \
              b * np.sin(theta_range) + \
              c * np.cos(2. * theta_range) + \
              d * np.sin(2. * theta_range)
        return res

    theta_range = np.arange(0., 2. * pi, theta_step)

    radius = radius_fct(theta_range, a, b, c, d)

    r_min = radius.min()

    radius += np.abs(r_min) + 0.2

    x_coords = radius * np.cos(theta_range)
    y_coords = radius * np.sin(theta_range)

    return x_coords, y_coords


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    x, y = random_closed_curve()

    plt.plot(x, y, 'r')
    plt.show()
