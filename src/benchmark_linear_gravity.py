"""
Run numerical experiments for quantum evolution
under gravitational potentials.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import matplotlib.lines as mlines
import time

from utils import create_grids
from initial_states import gaussian_wavepacket, normalize
from potentials import linear_gravity_potential
from solvers import time_evolution


# ------------------------------
# Global plotting style (CPC)
# ------------------------------

mpl.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 13,
    "legend.fontsize": 11,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "lines.linewidth": 2,
    "figure.figsize": (7, 5),
    "figure.dpi": 300,
    "axes.grid": False
})



# -----------------------------
# Physical parameters
# -----------------------------
hbar = 1.0
mass = 1.0
x0 = -10.0
p0 = 5.0
sigma = 1.0

#---------------------------------
# Number of branches
#----------------------------------
g1 = 0.8
g2 = 1.2
#g_list = [g1]      # single branch example
g_list = [g1, g2]     # multi-branch example   

B = len(g_list)


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

psi0 = np.zeros((B, N), dtype=complex)

for b in range(B):
     psi0[b] = gaussian_wavepacket(x, x0=x0, p0=p0, sigma=sigma, hbar=hbar)
     psi0[b] = normalize(psi0[b], dx)



# -----------------------------
# Potential
# -----------------------------

V_x = np.zeros((B, N))

for b, g in enumerate(g_list):
    V_x[b] = linear_gravity_potential(x, mass=mass, g=g)


start = time.perf_counter()

# -----------------------------
# Time evolution.   
# -----------------------------
result = time_evolution(
    psi0=psi0,
    V_x=V_x,
    x_grid=x,
    p_grid=p,
    dt=dt,
    n_steps=n_steps,
    hbar=hbar,
    mass=mass,
    store_wavefunction=True
)

end = time.perf_counter()
print("Runtime:", end - start)

# ---------------------------------
# Convergence tests
# ---------------------------------

# ============================================================
# Convergence test A: Linear Gravity (Physics Validation)
# ============================================================

# Compares each branch to same classical trajectory
# Returns global max error

def linear_gravity_error(dt_test, T_final,
                         psi0, x, p,
                         mass, hbar,
                         g1,
                         x0, p0,
                         time_evolution):

    B, N = psi0.shape

    n_steps = int(round(T_final / dt_test))

    # Build potential with correct shape
    V_x = np.zeros((B, N))
    for b in range(B):
        V_x[b] = mass * g1 * x

    result = time_evolution(
        psi0=psi0,
        V_x=V_x,
        x_grid=x,
        p_grid=p,
        dt=dt_test,
        n_steps=n_steps,
        hbar=hbar,
        mass=mass,
        store_wavefunction=False
    )

    t = result["time"]

    # shape (n_steps+1, B)
    x_num = result["x_expectation"]

    # Exact classical trajectory (scalar in time)
    x_exact = x0 + (p0 / mass) * t - 0.5 * g1 * t**2

    # Broadcast over branches
    error = np.max(np.abs(x_num - x_exact[:, None]))

    return error


# ============================================================
# Convergence test B: Nonlinear Self-Convergence (Order Validation)
# ============================================================

# Why max? Because, if one branch misbehaves, the scheme is not stable.

def nonlinear_self_convergence(dt_test, T_final,
                               psi0, x, p,
                               mass, hbar,
                               g1,
                               time_evolution):

    B, N = psi0.shape

    def run(dt):
        n_steps = int(round(T_final / dt))

        V_x = np.zeros((B, N))
        for b in range(B):
            V_x[b] = mass * g1 * x + 0.001 * x**3

        result = time_evolution(
            psi0=psi0,
            V_x=V_x,
            x_grid=x,
            p_grid=p,
            dt=dt,
            n_steps=n_steps,
            hbar=hbar,
            mass=mass,
            store_wavefunction=True
        )

        return result["psi"][-1]   # shape (B, N)

    psi_dt  = run(dt_test)
    psi_dt2 = run(dt_test / 2)
    psi_dt4 = run(dt_test / 4)

    # Compute error per branch
    E1 = np.linalg.norm(psi_dt  - psi_dt2, axis=1)
    E2 = np.linalg.norm(psi_dt2 - psi_dt4, axis=1)

    # Combine branches conservatively
    E1_max = np.max(E1)
    E2_max = np.max(E2)

    order = np.log2(E1_max / E2_max)

    return E1_max, order

# ============================================================
# Convergence tests: Main Execution
# ============================================================

if __name__ == "__main__":

    T_final = 5.0

    dt_list = [0.2, 0.1, 0.05, 0.025]

    print("\n--- Linear Gravity Test ---")
    linear_errors = []

    for dt in dt_list:
        err = linear_gravity_error(
            dt, T_final,
            psi0, x, p,
            mass, hbar,
            g1,
            x0, p0,
            time_evolution
        )
        print(f"dt = {dt:.5f}  error = {err:.3e}")
        linear_errors.append(err)

    print("\n--- Nonlinear Self-Convergence Test ---")
    dt_conv = []
    error_conv = []

    for dt in dt_list:
        err, order = nonlinear_self_convergence(
           dt, T_final,
           psi0, x, p,
           mass, hbar,
           g1,
           time_evolution
        )

        dt_conv.append(dt)
        error_conv.append(err)

        print(f"dt = {dt:.5f}   error = {err:.3e}    order = {order:.3e}")
        
        
# -------------------------------------------------
# Stability bound demonstration
# -------------------------------------------------
 
# compute t_max
g_stability = 8*g1
p_max = np.max(np.abs(p))
t_max = p_max / (mass * abs(g_stability))
print("Estimated spectral stability bound t_max =", t_max)

# define two simulations
T_safe = 0.36 * t_max
T_near = 0.95 * t_max
 
# error evolution function
def compute_error_vs_time(T_final, dt_test,
                          psi0, x, p,
                          mass, hbar,
                          g,
                          x0, p0,
                          time_evolution):

    B, N = psi0.shape
    n_steps = int(round(T_final / dt_test))

    # Build potential (B, N)
    V_x = np.zeros((B, N))
    for b in range(B):
        V_x[b] = mass * g * x

    result = time_evolution(
        psi0=psi0,
        V_x=V_x,
        x_grid=x,
        p_grid=p,
        dt=dt_test,
        n_steps=n_steps,
        hbar=hbar,
        mass=mass,
        store_wavefunction=False
    )

    t = result["time"]

    # shape (n_steps+1, B)
    x_num = result["x_expectation"]

    # Classical solution (time,)
    x_exact = x0 + (p0 / mass) * t - 0.5 * g * t**2

    # Broadcast over branches
    error_vals = np.abs(x_num - x_exact[:, None])

    # Conservative measure: worst branch
    error_max = np.max(error_vals, axis=1)

    return t, error_max
    
        
# Main execution
dt_test = 0.01

t_safe, err_safe = compute_error_vs_time(T_safe, dt_test, psi0, x, p, mass, hbar, g_stability, x0, p0, time_evolution)
t_near, err_near = compute_error_vs_time(T_near, dt_test, psi0, x, p, mass, hbar, g_stability, x0, p0, time_evolution)
err_safe[err_safe < 1e-16] = 1e-16
err_near[err_near < 1e-16] = 1e-16



#-------------------------------------------
# --- Figures placed in results file -----
#-------------------------------------------
output_dir = "results"
os.makedirs(output_dir, exist_ok=True)


# -------------------------------------------------
#Figure 1: CPC-style log norm deviation
# -------------------------------------------------
t = result["time"]
norm_dev = result["norm_deviation"]  # shape (T, B)

for b in range(norm_dev.shape[1]):
    plt.plot(t, norm_dev[:, b], label=f"Branch {b+1}")

plt.yscale('log')
plt.xlabel(r"$t$")
plt.ylabel(r'$|\|\psi_b(t)\|^2 - 1|$')
plt.grid(True, which="both", linestyle="--", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"{output_dir}/Gravity_potential_norm_deviation.png")
plt.close()

# -------------------------------------------------
# Figure 3: CPC-style log-log convergence plot
# -------------------------------------------------

dt_conv = np.array(dt_conv)
error_conv = np.array(error_conv)

# Sort for safety
idx = np.argsort(dt_conv)
dt_conv = dt_conv[idx]
error_conv = error_conv[idx]

# Reference slope-2 line (anchored to first point)
C_ref = error_conv[0] / dt_conv[0]**2
ref_line = C_ref * dt_conv**2


plt.loglog(dt_conv, error_conv,
           'o-', color='black',
           label=r'$\epsilon(\Delta t)$')

plt.loglog(dt_conv, ref_line,
           '--', color='black',
           label='Slope 2 reference')

plt.xlabel(r'$\Delta t$')
plt.ylabel(r'$\epsilon(\Delta t)$')

plt.grid(True, which="both", linestyle="--", alpha=0.3)
plt.legend()
plt.tight_layout()

plt.savefig(f"{output_dir}/nonlinear_convergence.png")
plt.close()


# --------------------------------------------
#  Figure 2: CPC-style Numerical vs Analytical Trajectory
# --------------------------------------------

# --- Plot trajectories ---

x_expect = result["x_expectation"]  # shape (T, B)

# Define branch colors (extendable if B > 2)
branch_colors = ["tab:blue", "tab:red", "tab:green", "tab:purple"]

for b in range(x_expect.shape[1]):

    color = branch_colors[b]

    # Numerical
    plt.plot(t, x_expect[:, b],
             color=color,
             linestyle="-")

    # Analytical (choose correct per branch)
    x_analytic = x0 + (p0 / mass) * t - 0.5 * g_list[b] * t**2

    plt.plot(t, x_analytic,
             color=color,
             linestyle="--")

plt.xlabel(r"$t$")
plt.ylabel(r"$\langle x(t) \rangle$")
plt.tick_params(labelsize=11)

g1_handle = mlines.Line2D([], [], color="tab:blue", linestyle="-",
                          label=r"$g_1$")
g2_handle = mlines.Line2D([], [], color="tab:red", linestyle="-",
                          label=r"$g_2$")

num_handle = mlines.Line2D([], [], color="black", linestyle="-",
                           label="Numerical")
ana_handle = mlines.Line2D([], [], color="black", linestyle="--",
                           label="Analytical")

plt.legend(handles=[g1_handle, g2_handle, num_handle, ana_handle],
           frameon=True,
           loc="best")
plt.grid(True, linestyle="--", alpha=0.3)
plt.tight_layout()
plt.savefig(f"{output_dir}/Linear_analytic_validation.png")
plt.close()

# -------------------------------------------------
# Figure 4: CPC-style Spectral Stability Bound Demonstration -- deviation vs time approaching t_max
# -------------------------------------------------


plt.plot(t_near, err_near, label="T = 0.95 t_max")
plt.plot(t_safe, err_safe, label="T = 0.36 t_max")
plt.axvline(t_max, linestyle="--", alpha=0.5)

plt.xlabel(r"$t$")
plt.ylabel(r"$|\langle x \rangle_{\mathrm{num}} - \langle x \rangle_{\mathrm{analytic}}|$")
plt.yscale("log")
plt.legend()
plt.grid(True, which="both", linestyle="--", alpha=0.3)
plt.tight_layout()
plt.savefig(f"{output_dir}/Gravity_potential_stability_bound.png")
plt.close()


# ---------------------------------------------------------
# Figure 5: CPC-style Branch norm conservation vs time
# ---------------------------------------------------------

# norm shape: (T, B)
norm = result["norm"]

# initial branch norms
norm0 = norm[0, :]   # shape (B,)

# deviation from initial value
dev = np.max(np.abs(norm - norm0), axis=1)

# Plot
plt.figure()
plt.plot(t, dev, label="max branch norm drift")
plt.yscale("log")
plt.xlabel(r"$t$")
plt.ylabel(r"$\max |N_b(t) - N_b(0)|$")
plt.grid(True, which="both", linestyle="--", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"{output_dir}/branch_norm_conservation.png")
plt.close()

# ---------------------------------------------------------
# Figure 6: CPC-style computational complexity
# ---------------------------------------------------------

# Grid sizes
N_val = np.array([1024, 2048, 4096, 8192, 16384])

# Measured runtimes (seconds)
runtime_B1 = np.array([0.16, 0.24, 0.39, 0.72, 1.35])      # For B=1
runtime_B2 = np.array([0.26, 0.38, 0.68, 1.27, 2.69])      # For B=2

# Compute N log2 N
N_logN = N_val * np.log2(N_val)

# Scaling diagnostic
scaled_runtime_1 = runtime_B1 / N_logN
scaled_runtime_2 = runtime_B2 / N_logN


# Create figure with two panels
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10,4))

# Left plot: runtime vs N log N
ax1.plot(N_logN, runtime_B1, 'o-', linewidth=2, markersize=6, label='B = 1')
ax1.plot(N_logN, runtime_B2, 's-', linewidth=2, markersize=6, label='B = 2')

ax1.set_xlabel(r"$N \log_2 N$ (grid scaling variable)", fontsize=12)
ax1.set_ylabel("Runtime (s)", fontsize=12)
ax1.grid(True, linestyle='--', alpha=0.3)
ax1.legend()

# Right plot: normalized runtime
ax2.plot(N_val, scaled_runtime_1, 'o-', linewidth=2, markersize=6, label='B = 1')
ax2.plot(N_val, scaled_runtime_2, 's-', linewidth=2, markersize=6, label='B = 2')

ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
ax2.yaxis.get_offset_text().set_visible(False)
ax2.set_ylabel(r"$T/(N\log_2 N) \times 10^{-5}$")

ax2.set_xlabel(r"$N$", fontsize=12)
#ax2.set_ylabel(r"$T/(N\log_2 N)$", fontsize=12)
ax2.grid(True, linestyle='--', alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.savefig(f"{output_dir}/comp_complexity.png")
plt.close()
