"""
Observable calculations.
"""

import numpy as np

def expectation_x(psi, x):
    """
    Expectation value ⟨x⟩.
    """
    dx = x[1] - x[0]
    return np.sum(np.abs(psi)**2 * x) * dx

# Operator-based expectation

def expectation_p_op(psi, x, hbar):
    """
    Expectation value ⟨p⟩ computed in position space
    using the momentum operator.
    """
    dx = x[1] - x[0]

    # Central finite difference for derivative
    dpsi_dx = np.gradient(psi, dx)

    return np.real(
        np.sum(np.conj(psi) * (-1j * hbar * dpsi_dx)) * dx
    )


# FFT momentum expectation

def expectation_p_fft(psi, p):
    """
    Expectation value ⟨p⟩ using FFT (with correct measure).
    """
    psi_p = np.fft.fft(psi)

# fftshift is required to align momentum grid with Fourier components

    psi_p = np.fft.fftshift(psi_p)
    p_shifted = np.fft.fftshift(p)

    dp = p_shifted[1] - p_shifted[0]

    prob_density_p = np.abs(psi_p)**2

    # Normalize with measure
    norm = np.sum(prob_density_p) * dp
    prob_density_p /= norm

    return np.sum(p_shifted * prob_density_p) * dp



def overlap(psi1, psi2, dx):
    """
    Overlap |<psi1 | psi2>|.
    """
    raise NotImplementedError
