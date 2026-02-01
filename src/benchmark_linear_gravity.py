"""
Run numerical experiments for quantum evolution
under gravitational potentials.
"""

import numpy as np
import matplotlib.pyplot as plt

from utils import create_grids
from initial_states import gaussian_wavepacket, normalize
from potentials import linear_gravity_potential
from solvers import time_evolution


# -----------------------------
# Physical parameters
# -----------------------------
hbar = 1.0
mass = 1.0
x0 = -10.0
p0 = 5.0
sigma = 1.0

g1 = 0.8
g2 = 1.2
geff = 1

# -----------------------------
# Numerical parameters
# -----------------------------
L = 100.0
N = 2048

dt = 0.01
n_steps = 1000

# -----------------------------
# Create grids
# -----------------------------
x, p, dx = create_grids(N, L, hbar)

# -----------------------------
# Initial state
# -----------------------------
psi0 = gaussian_wavepacket(x, x0=x0, p0=p0, sigma=sigma, hbar=hbar)
psi0 = normalize(psi0, dx)


# -----------------------------
# Potential
# -----------------------------
V1 = linear_gravity_potential(x, mass=mass, g=g1)
V2 = linear_gravity_potential(x, mass=mass, g=g2)
Veff = linear_gravity_potential(x, mass=mass, g=geff)


# -----------------------------
# Time evolution
# -----------------------------
results_1 = time_evolution(
    psi0=psi0,
    V_x=V1,
    x_grid=x,
    p_grid=p,
    dt=dt,
    n_steps=n_steps,
    hbar=hbar,
    mass=mass,
    store_wavefunction=True
)

results_2 = time_evolution(
    psi0=psi0,
    V_x=V2,
    x_grid=x,
    p_grid=p,
    dt=dt,
    n_steps=n_steps,
    hbar=hbar,
    mass=mass,
    store_wavefunction=True
)

results_eff = time_evolution(
    psi0=psi0,
    V_x=Veff,
    x_grid=x,
    p_grid=p,
    dt=dt,
    n_steps=n_steps,
    hbar=hbar,
    mass=mass,
    store_wavefunction=False
)

# -----------------------------
# Analytic solutions
# -----------------------------
t = results_1["time"]
x1_analytic = x0 + (p0 / mass) * t - 0.5 * g1 * t**2  
p1_analytic = p0 - mass * g1 * t  
x2_analytic = x0 + (p0 / mass) * t - 0.5 * g2 * t**2  
p2_analytic = p0 - mass * g2 * t  
xeff_analytic = x0 + (p0 / mass) * t - 0.5 * geff * t**2  
peff_analytic = p0 - mass * geff * t 


# =====================================================
# Overlap (Case A only)
# =====================================================
psi1_t = results_1["psi"]
psi2_t = results_2["psi"]

overlap = np.array([
    np.abs(np.trapz(np.conj(psi1_t[n]) * psi2_t[n], x))
    for n in range(n_steps)
])


# -----------------------------
# Plots
# -----------------------------
plt.figure()
plt.plot(t, results_1["x_expectation"], label="Numerical ⟨x⟩_1")
plt.plot(t, results_2["x_expectation"], label="Numerical ⟨x⟩_2")
plt.plot(t, results_eff["x_expectation"], label="Numerical ⟨x⟩_eff")
plt.plot(t, x1_analytic, "--", label="Analytic ⟨x⟩_1")
plt.plot(t, x2_analytic, "--", label="Analytic ⟨x⟩_2")
plt.plot(t, xeff_analytic, "--", label="Analytic ⟨x⟩_eff")
plt.xlabel("Time")
plt.ylabel("⟨x⟩")
plt.legend()
plt.title("Gravity Potential: Position Expectation")
plt.show()

plt.figure()
plt.plot(t, results_1["p_expectation"], label="Numerical ⟨p⟩_1")
plt.plot(t, results_2["p_expectation"], label="Numerical ⟨p⟩_2")
plt.plot(t, results_eff["p_expectation"], label="Numerical ⟨p⟩_eff")
plt.plot(t, p1_analytic, "--", label="Analytic ⟨p⟩_1")
plt.plot(t, p2_analytic, "--", label="Analytic ⟨p⟩_2")
plt.plot(t, peff_analytic, "--", label="Analytic ⟨p⟩_eff")
plt.xlabel("Time")
plt.ylabel("⟨p⟩")
plt.legend()
plt.title("Gravity Potential: Momentum Expectation")
plt.show()

plt.figure()
plt.plot(t, overlap)
plt.xlabel("t")
plt.ylabel(r"$|\langle \psi_1 | \psi_2 \rangle|$")
plt.title("Branch overlap")
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(t, results_1["norm"], label="Branch 1")
plt.plot(t, results_2["norm"], label="Branch 2")
plt.plot(t, results_eff["norm"], label="Effective")
plt.xlabel("Time")
plt.ylabel("Norm")
plt.title("Wavefunction Norm")
plt.show()
