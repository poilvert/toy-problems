#!/usr/bin/env python

import matplotlib.pyplot as plt
from curvature_flow_1d import curve_flow
from generate_random_closed_curve import random_closed_curve

try:
    curve_flow = profile(curve_flow)
except NameError:
    pass

# -- generate random closed curve in the cartesian plane
x, y = random_closed_curve(a=0.12, b=-0.25, c=-0.41, d=0.37)

# -- evolve the shape of the random closed curve with
#    curvature flow
x_l, y_l = curve_flow(x, y, dt=0.01, nstep=400)

# -- plot evolution of shape
start = 0
fourth = len(x_l) / 4
half = len(x_l) / 2
three_fourth = 3 * len(x_l) / 4
end = len(x_l) - 1
plt.plot(x_l[start], y_l[start], '-',
         x_l[fourth], y_l[fourth],'--',
         x_l[half], y_l[half],'-.',
         x_l[three_fourth], y_l[three_fourth], '.',
         x_l[end], y_l[end], '*')
plt.show()
