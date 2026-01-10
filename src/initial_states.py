"""
Initial quantum states.
"""

import numpy as np

def gaussian_wavepacket(x, x0, p0, sigma, hbar):
    """
    Gaussian wave packet in position space.

    Parameters
    ----------
    x : ndarray
        Spatial grid.
    x0 : float
        Initial position expectation value.
    p0 : float
        Initial momentum expectation value.
    sigma : float
        Width of the packet.
    hbar : float
        Reduced Planck constant.

    Returns
    -------
    psi0 : ndarray (complex)
        Normalized wavefunction ψ(x,0).
    """
    psi0 = np.exp(- (x - x0)**2 / (4.0 * sigma**2)+ 1j * p0 * (x - x0) / hbar)
    
    """
    Normalization of Gaussian in the continuum - not used
    """
    # normalization = (1.0 / (2.0 * np.pi * sigma**2))**0.25
    # psi0 = normalization * psi0    

    return psi0



"""
Numerical normalization of Gaussian - used
"""

def normalize(psi, dx):
    """
    Normalize wavefunction on a discrete grid.
    """
    return psi / np.sqrt(np.sum(np.abs(psi)**2) * dx)



