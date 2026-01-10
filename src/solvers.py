"""
Numerical solvers for the time-dependent Schrödinger equation.
"""

import numpy as np

def split_operator_step(psi_x, V_x, U_T, dt, hbar):
    """
    Perform a single split-operator time step.

    Parameters
    ----------
    psi_x : ndarray (complex)
        Wavefunction in position space at time t.
    V_x : ndarray
        Potential evaluated on the spatial grid.
    U_T : ndarray
        Precomputed kinetic evolution operator in momentum space.
    dt : float
        Time step.
    hbar : float

    Returns
    -------
    psi_x_new : ndarray (complex)
        Wavefunction after one time step.
    """

    # --- Half-step potential evolution (position space)
    U_V_half = np.exp(-1j * V_x * dt / (2.0 * hbar))
    psi_x = U_V_half * psi_x

    # --- Transform to momentum space
    psi_p = np.fft.fft(psi_x)

    # --- Full-step kinetic evolution (momentum space)
    psi_p = U_T * psi_p

    # --- Transform back to position space
    psi_x = np.fft.ifft(psi_p)

    # --- Half-step potential evolution (position space)
    psi_x = U_V_half * psi_x

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
    Time evolution using the split-operator Fourier method.

    Parameters
    ----------
    psi0 : ndarray (complex)
        Initial wavefunction in position space.
    V_x : ndarray
        Potential evaluated on the spatial grid.
    x_grid : ndarray
        Spatial grid.
    p_grid : ndarray
        Momentum grid.
    dt : float
        Time step.
    n_steps : int
        Number of time steps.
    hbar : float
        Reduced Planck constant.
    mass : float, optional
        Particle mass.
    store_wavefunction : bool, optional
        If True, store full wavefunction at each time step.

    Returns
    -------
    results : dict
        Dictionary containing time series of observables
        and optionally wavefunctions.
    """

    # --- Grid spacing
    dx = x_grid[1] - x_grid[0]

    # --- Precompute kinetic evolution operator
    U_T = np.exp(-1j * p_grid**2 * dt / (2.0 * mass * hbar))

    # --- Initialize storage
    times = np.arange(n_steps) * dt
    x_expect = np.zeros(n_steps)
    p_expect = np.zeros(n_steps)
    norm = np.zeros(n_steps)

    if store_wavefunction:
        psi_t = np.zeros((n_steps, len(psi0)), dtype=complex)

    # --- Initialize wavefunction
    psi = psi0.copy()

    # --- Time evolution loop
    for n in range(n_steps):

        # Store observables
        x_expect[n] = expectation_x(psi, x_grid)
        p_expect[n] = expectation_p_op(psi, x_grid, hbar)
        norm[n] = np.sum(np.abs(psi)**2) * dx
        

        if store_wavefunction:
            psi_t[n] = psi

        # Propagate one time step
        psi = split_operator_step(psi, V_x, U_T, dt, hbar)

    # --- Collect results
    results = {
        "time": times,
        "x_expectation": x_expect,
        "p_expectation": p_expect,
        "norm": norm
    }

    if store_wavefunction:
        results["psi"] = psi_t

    return results
