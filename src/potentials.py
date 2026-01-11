"""
Gravitational potentials.
"""

import numpy as np

def free_particle(x):
    """
    Free particle potential V(x) = 0.

    Parameters
    ----------
    x : ndarray
        Spatial grid.

    Returns
    -------
    V : ndarray
        Zero potential evaluated on the grid.
    """
    return np.zeros_like(x)


def linear_gravity_potential(x, mass, g):
    """
    Linear gravitational potential V(x) = m g x.
    """
    return mass * g * x
