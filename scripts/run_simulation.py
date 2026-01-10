"""
Run numerical experiments for quantum evolution
under gravitational potentials.
"""

import numpy as np
from src.solvers import time_evolution
from src.potentials import linear_gravitational_potential
from src.initial_states import gaussian_wavepacket
from src.utils import spatial_grid, momentum_grid

# Physical parameters
hbar = 1.0
mass = 1.0
g = 1.0

# Numerical parameters
L = 50.0
N = 1024
dt = 0.01
n_steps = 1000

# Grids
x, dx = spatial_grid(L, N)
p = momentum_grid(N, dx, hbar)

# Initial state
psi0 = gaussian_wavepacket(
    x=x,
    x0=0.0,
    sigma=1.0,
    p0=0.0,
    hbar=hbar
)

# Potential
V_x = linear_gravitational_potential(x, mass, g)

# Time evolution
results = time_evolution(
    psi0=psi0,
    V_x=V_x,
    x_grid=x,
    p_grid=p,
    dt=dt,
    n_steps=n_steps,
    hbar=hbar
)
