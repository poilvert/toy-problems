#!/usr/bin/env python

__author__ = 'Nicolas Poilvert'

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
from scipy.spatial.distance import cdist

# -----------------------------------------------
# Main Molecular Dynamics (Lennard-Jones) routine
# -----------------------------------------------
def molecular_dynamics(
        nstep=5000,
        dt=0.5,
        m=39.94,
        sigma=3.4,
        epsilon=0.0104):

    # -- initializing positions
    xs, ys = get_initial_positions()
    assert xs.ndim == ys.ndim
    assert xs.ndim == 1
    assert xs.size > 1
    natoms = xs.size

    # -- basic tensors for the MD run
    ndim = 2
    r = np.zeros((natoms, ndim, nstep+2), dtype='d')
    v = np.zeros((natoms, ndim, nstep+2), dtype='d')

    # -- initializing positions
    r[:, :, 0] = np.array([xs, ys]).T
    r[:, :, 1] = r[:, :, 0]

    # -- loop invariants (useful constants used during the MD run)
    acc_prefactor = dt ** 2 / m
    vel_prefactor = 1. / (2. * dt)

    # -- visualizing the particles
    plt.ion()
    plt.xlabel('x coordinate', fontsize=16)
    plt.ylabel('y coordinate', fontsize=16)
    plt.plot(r[:, 0, 0], r[:, 1, 0], 'bo')
    plt.title('time step %04i / %04i' % (0, nstep))

    # -- main loop (MD run itself)
    T, U, KE = [], [], []
    for t in range(1, nstep+1):

        # -- Computing the Lennard-Jones forces
        f = LJ_forces(r[..., t], sigma, epsilon)

        # -- Computing the Lennard-Jones energy
        U += [LJ_energy(r[..., t], sigma, epsilon)]

        # -- Verlet Algorithm
        r[..., t+1] = 2. * r[..., t] - r[..., t-1] + acc_prefactor * f
        v[..., t] = vel_prefactor * (r[..., t+1] - r[..., t-1])

        # -- Computing the Kinetic Energy
        KE += [Kinetic_energy(v[..., t], m)]
        T += [0.0102 * t*dt]

        # -- Update visualization
        if t%25 == 0:
            plt.plot(r[:, 0, t], r[:, 1, t], 'bo')
            plt.title('time step %04i / %04i' % (t, nstep))
            plt.draw()
            plt.pause(0.05)

    plt.clf()

    return np.array(T), np.array(U), np.array(KE)

# --------------------------------------------
# Utility routines used by the main MD routine
# --------------------------------------------
def LJ_forces(r, sigma, epsilon):
    natoms = r.shape[0]
    pairwise_distances = cdist(r, r)
    pairwise_distances[np.diag_indices(natoms)] = 1.
    s_over_r6 = (sigma / pairwise_distances) ** 6
    pairwise_scalars = 4.*epsilon * (12.*s_over_r6 ** 2 - 6.*s_over_r6)
    pairwise_scalars *= 1. / pairwise_distances ** 2
    pairwise_vectors = r[:, None, :] - r[None, :, :]
    return np.sum(pairwise_scalars[..., None] * pairwise_vectors, axis=1)

def LJ_energy(r, sigma, epsilon):
    pairwise_distances = pdist(r)
    s_over_r6 = (sigma / pairwise_distances) ** 6
    return 4*epsilon * np.sum(s_over_r6 ** 2 - s_over_r6)

def Kinetic_energy(v, m):
    return 0.5 * m * np.sum(v**2)

# ----------------------------------------------------------
# Routines used to select the initial positions of the atoms
# ----------------------------------------------------------
class AtomBuilder:
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        if event.inaxes != self.line.axes:
            return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        self.line.set_linestyle('None')
        self.line.set_marker('o')
        self.line.set_markerfacecolor('b')
        self.line.figure.canvas.draw()

    def dump_coordinates(self):
        return np.array(self.xs), np.array(self.ys)

def get_initial_positions(R=20.):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title('Click to choose atomic positions, then close window')
    ax.set_xlim((-R, R))
    ax.set_ylim((-R, R))
    ax.set_xlabel('x coordinate')
    ax.set_ylabel('y coordinate')
    line, = ax.plot([], [])
    builder = AtomBuilder(line)
    plt.show()
    return builder.dump_coordinates()


if __name__ == "__main__":

    # -- Performing the MD run
    t, U, KE = molecular_dynamics()

    # -- Visualization of the energy terms along the MD path
    plt.ioff()
    plt.plot(t, U, 'ro-', linewidth=2, label='potential energy U')
    plt.plot(t, KE, 'b*-', linewidth=2, label='kinetic energy K')
    plt.plot(t, U+KE, 'k^-', linewidth=2, label='total energy E')
    plt.xlabel('time t (ps)', fontsize=16)
    plt.ylabel('energy (eV)', fontsize=16)
    plt.xlim((0., t.max()))
    plt.legend(loc='best')
    plt.show()
