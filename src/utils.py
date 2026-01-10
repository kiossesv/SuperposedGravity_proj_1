"""
Utility functions.
"""

import numpy as np

def create_grids(N, L, hbar=1.0):
    """
    Create consistent position and momentum grids.
    """
    x = np.linspace(-L/2, L/2, N, endpoint=False)
    dx = x[1] - x[0]

    p = np.fft.fftfreq(N, d=dx) * 2 * np.pi * hbar

    return x, p, dx
