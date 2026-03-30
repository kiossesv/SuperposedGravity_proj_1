"""
Observable calculations.
"""

import numpy as np

def expectation_x(psi, x):
    """
    Expectation value ⟨x⟩ per branch.
    psi shape: (B, N)
    Returns: shape (B,)
    """
    assert psi.ndim == 2, "psi must have shape (B, N)"

    dx = x[1] - x[0]

    prob_density = np.abs(psi)**2

    return np.sum(prob_density * x, axis=1) * dx
    
    
    
# Operator-based momentum expectation

def expectation_p_op(psi, x, hbar):
    """
    Expectation value ⟨p⟩ per branch using
    momentum operator in position space.
    psi shape: (B, N)
    Returns: shape (B,)
    """
    assert psi.ndim == 2, "psi must have shape (B, N)"

    dx = x[1] - x[0]

    # Derivative along spatial dimension only
    dpsi_dx = np.gradient(psi, dx, axis=1)

    return np.real(
        np.sum(np.conj(psi) * (-1j * hbar * dpsi_dx), axis=1) * dx
    )
    
    


def overlap(psi1, psi2, dx):
    """
    Overlap |<psi1 | psi2>|.
    """
    raise NotImplementedError
