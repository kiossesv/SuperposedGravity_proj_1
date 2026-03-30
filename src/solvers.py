"""
Numerical solvers for the time-dependent Schrödinger equation.
"""

import numpy as np

def split_operator_step(psi_x, V_x, U_T, dt, hbar):
    """
    Perform a single split-operator time step
    for B independent branches.

    Parameters
    ----------
    psi_x : ndarray (complex), shape (B, N)
    V_x   : ndarray, shape (B, N)
    U_T   : ndarray, shape (N,)
    dt    : float
    hbar  : float

    Returns
    -------
    psi_x_new : ndarray (complex), shape (B, N)
    """

    # --- Sanity checks (good for CPC-level robustness)
    assert psi_x.ndim == 2, "psi_x must have shape (B, N)"
    assert V_x.shape == psi_x.shape, "V_x must match psi_x shape"
    assert U_T.ndim == 1, "U_T must have shape (N,)"

    # --- Half-step potential
    U_V_half = np.exp(-1j * V_x * dt / (2.0 * hbar))
    psi_x = U_V_half * psi_x

    # --- FFT along spatial dimension only
    psi_p = np.fft.fft(psi_x, axis=1)

    # --- Kinetic evolution
    psi_p *= U_T

    # --- Inverse FFT
    psi_x = np.fft.ifft(psi_p, axis=1)

    # --- Half-step potential
    psi_x *= U_V_half

    return psi_x



from observables import expectation_x, expectation_p_op

def time_evolution(
    psi0,
    V_x,
    x_grid,
    p_grid,
    dt,
    n_steps,
    hbar,
    mass,
    store_wavefunction=False
):
    """
    Time evolution using split-operator method
    for B independent branches.
    """

    # --- Enforce (B, N)
    assert psi0.ndim == 2, "psi0 must have shape (B, N)"
    assert V_x.shape == psi0.shape, "V_x must match psi0 shape"

    B, N = psi0.shape

    dx = x_grid[1] - x_grid[0]

    # --- Precompute kinetic operator (1D)
    U_T = np.exp(-1j * p_grid**2 * dt / (2.0 * mass * hbar))

    # --- Storage arrays (per branch)
    times = np.linspace(0, n_steps * dt, n_steps + 1)

    x_expect = np.zeros((n_steps + 1, B))
    p_expect = np.zeros((n_steps + 1, B))
    norm = np.zeros((n_steps + 1, B))

    if store_wavefunction:
        psi_t = np.zeros((n_steps + 1, B, N), dtype=complex)

    # --- Initialize wavefunction
    psi = psi0.copy()

    # --- Time evolution loop
    for n in range(n_steps):

        # Observables per branch
        x_expect[n] = expectation_x(psi, x_grid)
        p_expect[n] = expectation_p_op(psi, x_grid, hbar)

        norm[n] = np.sum(np.abs(psi)**2, axis=1) * dx

        if store_wavefunction:
            psi_t[n] = psi

        # Propagate
        psi = split_operator_step(psi, V_x, U_T, dt, hbar)

    # --- Final state storage
    x_expect[-1] = expectation_x(psi, x_grid)
    p_expect[-1] = expectation_p_op(psi, x_grid, hbar)
    norm[-1] = np.sum(np.abs(psi)**2, axis=1) * dx

    if store_wavefunction:
        psi_t[-1] = psi

    # --- Norm deviation per branch
    norm_deviation = np.abs(norm - 1.0)
    norm_deviation[norm_deviation < 1e-16] = 1e-16

    # --- Collect results
    results = {
        "time": times,
        "x_expectation": x_expect,
        "p_expectation": p_expect,
        "norm": norm,
        "norm_deviation": norm_deviation
    }

    if store_wavefunction:
        results["psi"] = psi_t

    return results