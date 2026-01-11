import numpy as np
import matplotlib.pyplot as plt

from utils import create_grids
from initial_states import gaussian_wavepacket, normalize
from potentials import free_particle
from solvers import time_evolution

# -----------------------------
# Physical parameters
# -----------------------------
hbar = 1.0
mass = 1.0

x0 = -10.0
p0 = 2.0
sigma = 1.0

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
V = free_particle(x)

# -----------------------------
# Time evolution
# -----------------------------
results = time_evolution(
    psi0=psi0,
    V_x=V,
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
t = results["time"]
x_analytic = x0 + (p0 / mass) * t
p_analytic = p0 * np.ones_like(t)

# -----------------------------
# Plots
# -----------------------------
plt.figure()
plt.plot(t, results["x_expectation"], label="Numerical ⟨x⟩")
plt.plot(t, x_analytic, "--", label="Analytic ⟨x⟩")
plt.xlabel("Time")
plt.ylabel("⟨x⟩")
plt.legend()
plt.title("Free Particle: Position Expectation")
plt.show()

plt.figure()
plt.plot(t, results["p_expectation"], label="Numerical ⟨p⟩")
plt.plot(t, p_analytic, "--", label="Analytic ⟨p⟩")
plt.xlabel("Time")
plt.ylabel("⟨p⟩")
plt.legend()
plt.title("Free Particle: Momentum Expectation")
plt.show()


plt.figure()
plt.plot(t, results["norm"])
plt.xlabel("Time")
plt.ylabel("Norm")
plt.title("Wavefunction Norm")
plt.show()
