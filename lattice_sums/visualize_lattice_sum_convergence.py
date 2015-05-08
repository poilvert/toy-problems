#!/usr/bin/env python

# Authors: Nicolas Poilvert
# License: BSD

"""
A way to visualize the efficiency of the Ewald summation technique compared to
direct summations when it comes to computing lattice sums.
"""

import numpy as np
from cubic_3D_lattices import lattice_sum, Ewald_lattice_sum
import matplotlib.pyplot as plt


def visualize_convergence(l, m, crystal_type='sc'):
    # -- increasing size of the lattice grid
    N_grid = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 30]
    N_grid_points = [(2*i+1)**3 - 1 for i in N_grid]
    # -- direct summation
    direct_sums = [
            lattice_sum(l, m, crystal_type=crystal_type,
                N=N)
            for N in N_grid
        ]
    Ewald_sums = [
            Ewald_lattice_sum(l, m, crystal_type=crystal_type,
                N=N)
            for N in N_grid
        ]
    # -- We take the Ewald sum value for the largest N as a reference for what
    # the lattice sum value should be when converged
    ref_value = Ewald_sums[-1]
    direct_errors = np.abs(ref_value - np.array(direct_sums))
    plt.plot(N_grid_points, direct_errors,
             'ro-', label='Direct', linewidth=3)
    Ewald_errors = np.abs(ref_value - np.array(Ewald_sums))
    plt.plot(N_grid_points, Ewald_errors,
             'bo-', label='Ewald', linewidth=3)
    plt.xscale('log')
    plt.xlim((0., max(N_grid_points)))
    plt.ylim((0., max(direct_errors.max(), Ewald_errors.max())))
    plt.xlabel('# of grid points for summation', fontsize=16)
    plt.ylabel('Absolute value of error', fontsize=16)
    plt.legend(loc='best')
    plt.title('Crystal Type : %s, (l, m) : (%i, %i)' % (crystal_type.upper(),
        l, m))
    plt.show()

if __name__ == "__main__":

    # -- type of cubic crystal to test ('sc', 'bcc' or 'fcc')
    crystal_type = 'sc'

    # -- order of the lattice sum (see the definition of the lattice sum to know
    # what those parameters mean)
    l, m = 4, 0

    # -- plotting both the direct summation and Ewald summation to compare them
    visualize_convergence(l, m, crystal_type)
