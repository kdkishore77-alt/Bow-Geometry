"""
Reproduces all figures and sensitivity analyses for the manuscript.

Run:
    python src/main.py
"""

import numpy as np
import matplotlib.pyplot as plt

# Deterministic behavior
np.random.seed(0)

from parameters import *
from geometry import *
from forces import *
from utils import *
from simulations import *


# =====================================
# MAXIMUM DRAW ANGLES
# =====================================

theta_max_self = find_theta_max(
    ell_self, y0_self, L0_self, theta0_self, X_MAX_SELF
)
theta_max_rec = find_theta_max(
    ell_rec, y0_rec, L0_rec, theta0_rec, X_MAX_REC
)
theta_max_yumi = find_theta_max(
    ell_avg_yumi, y0_yumi, L0_yumi, theta0_yumi, X_MAX_YUMI
)


# =====================================
# ANGLE GRIDS
# =====================================

theta_self = np.linspace(theta0_self + 1e-4, theta_max_self, 600)
theta_rec  = np.linspace(theta0_rec  + 1e-4, theta_max_rec,  600)
theta1_yumi = np.linspace(theta0_yumi + 1e-4, theta_max_yumi, 600)

theta2_yumi = solve_yumi_secondary_angle(theta1_yumi)

# =====================================
# DRAW LENGTHS
# =====================================

x_self = draw_length(theta_self, ell_self, y0_self, L0_self, theta0_self)
x_rec  = draw_length(theta_rec,  ell_rec,  y0_rec,  L0_rec,  theta0_rec)

x_yumi = (
    draw_length(theta1_yumi, ell1_yumi, y0_yumi, L0_yumi, theta0_yumi)
  + draw_length(theta2_yumi, ell2_yumi, y0_yumi, L0_yumi, theta0_yumi)
)

# =====================================
# STRING ANGLES
# =====================================

ph_self = phi(theta_self, ell_self, y0_self, L0_self)
ph_rec  = phi(theta_rec,  ell_rec,  y0_rec,  L0_rec)
ph1 = phi(theta1_yumi, ell1_yumi, y0_yumi, L0_yumi)
ph2 = phi(theta2_yumi, ell2_yumi, y0_yumi, L0_yumi)

dx_dtheta_self = np.gradient(x_self, theta_self)
dx_dtheta_rec  = np.gradient(x_rec,  theta_rec)
dx_dtheta_yumi = np.gradient(x_yumi, theta1_yumi)


# =====================================
# GEOMETRIC ADMISSIBILITY
# =====================================

mask_self = (np.abs(np.cos(ph_self - theta_self)) > COS_EPS) & (dx_dtheta_self > 0)
mask_rec  = (np.abs(np.cos(ph_rec  - theta_rec )) > COS_EPS) & (dx_dtheta_rec  > 0)
mask_yumi = (np.abs(np.cos(ph1 - theta1_yumi)) > COS_EPS) & (dx_dtheta_yumi > 0)

theta_self, x_self = theta_self[mask_self], x_self[mask_self]
theta_rec,  x_rec  = theta_rec[mask_rec],  x_rec[mask_rec]

theta1_yumi = theta1_yumi[mask_yumi]
theta2_yumi = theta2_yumi[mask_yumi]
x_yumi      = x_yumi[mask_yumi]


# =====================================
# FORCES
# =====================================

F_self_x = F_self(theta_self)
F_rec_x  = F_recurve(theta_rec)

ph1 = phi(theta1_yumi, ell1_yumi, y0_yumi, L0_yumi)
ph2 = phi(theta2_yumi, ell2_yumi, y0_yumi, L0_yumi)

rc1 = np.maximum(ell1_yumi * np.cos(ph1 - theta1_yumi), 1e-4)
rc2 = np.maximum(ell2_yumi * np.cos(ph2 - theta2_yumi), 1e-4)

T_yumi = k1_yumi * (theta1_yumi - theta0_yumi) / rc1
F_yumi_x = T_yumi * (np.cos(ph1) + np.cos(ph2))


# =====================================
# FIGURE SETUP
# =====================================

fig, axes = plt.subplots(3, 2, figsize=(8.27, 11.69), constrained_layout=True)

# =====================================
# (a) FORCE–DRAW CURVES
# =====================================

ax = axes[0, 0]
ax.plot(x_self, F_self_x, lw=3, label="Self-Bow")
ax.plot(x_rec,  F_rec_x,  lw=3, label="Composite Recurve")
ax.plot(x_yumi, F_yumi_x, lw=3, label="Asymmetric Yumi")
ax.set_xlabel(r"Draw length $x$ (m)", fontsize=16)
ax.set_ylabel(r"Draw force $F$ (N)", fontsize=16)
ax.legend()
ax.set_title("(a)", loc="left", fontweight="bold")

# =====================================
# (b) STORED ENERGY
# =====================================

E_self = np.array([stored_energy(F_self_x[:i], x_self[:i]) for i in range(10, len(x_self))])
E_rec  = np.array([stored_energy(F_rec_x[:i],  x_rec[:i])  for i in range(10, len(x_rec))])
E_yumi = np.array([stored_energy(F_yumi_x[:i], x_yumi[:i]) for i in range(10, len(x_yumi))])

ax = axes[0, 1]
ax.plot(x_self[10:], E_self, lw=3, label="Self-Bow")
ax.plot(x_rec[10:],  E_rec,  lw=3, label="Composite Recurve")
ax.plot(x_yumi[10:], E_yumi, lw=3, label="Asymmetric Yumi")
ax.set_xlabel("Draw length $x$ (m)", fontsize=16)
ax.set_ylabel("Stored energy $U$ (J)", fontsize=16)
ax.legend()
ax.set_title("(b)", loc="left", fontweight="bold")

# =====================================
# (c) STATIC EFFICIENCY VS BRACE ANGLE
# =====================================

theta0_vals = np.linspace(10*np.pi/180, 40*np.pi/180, 50)
eta_self, eta_rec, eta_yumi = [], [], []

for th0 in theta0_vals:

    L0_temp = y0_self + ell_self * np.cos(th0)
    theta_max_temp = find_theta_max(
        ell_self, y0_self, L0_temp, th0, X_MAX_SELF
    )
    theta_temp = np.linspace(th0 + 1e-4, theta_max_temp, 400)
    ph_temp = phi(theta_temp, ell_self, y0_self, L0_temp)
    x_temp = draw_length(theta_temp, ell_self, y0_self, L0_temp, th0)
    F_s = (2 * k_self * (theta_temp - th0) * np.cos(ph_temp)) / \
          (ell_self * np.cos(ph_temp - theta_temp))
    eta_self.append(static_efficiency(F_s, x_temp))

    L0_temp = y0_rec + ell_rec * np.cos(th0)
    theta_max_temp = find_theta_max(
        ell_rec, y0_rec, L0_temp, th0, X_MAX_REC
    )
    theta_temp = np.linspace(th0 + 1e-4, theta_max_temp, 400)
    ph_temp = phi(theta_temp, ell_rec, y0_rec, L0_temp)
    x_temp = draw_length(theta_temp, ell_rec, y0_rec, L0_temp, th0)
    F_r = (2 * k_rec * (theta_temp - th0) * np.cos(ph_temp)) / \
          (ell_rec * np.cos(ph_temp - theta_temp) +
           ell_r * np.exp(-alpha * (theta_temp - th0)))
    eta_rec.append(static_efficiency(F_r, x_temp))

    ell_avg = ell_avg_yumi
    k_avg   = k_avg_yumi
    L0_temp = y0_yumi + ell_avg * np.cos(th0)
    theta_max_temp = find_theta_max(
        ell_avg, y0_yumi, L0_temp, th0, X_MAX_YUMI
    )
    theta_temp = np.linspace(th0 + 1e-4, theta_max_temp, 400)
    ph_temp = phi(theta_temp, ell_avg, y0_yumi, L0_temp)
    x_temp = draw_length(theta_temp, ell_avg, y0_yumi, L0_temp, th0)
    F_y = (2 * k_avg * (theta_temp - th0) * np.cos(ph_temp)) / \
          (ell_avg * np.cos(ph_temp - theta_temp))
    eta_yumi.append(static_efficiency(F_y, x_temp))

ax = axes[1, 0]
ax.plot(theta0_vals*180/np.pi, eta_self, lw=3, label="Self-Bow")
ax.plot(theta0_vals*180/np.pi, eta_rec,  lw=3, label="Composite Recurve")
ax.plot(theta0_vals*180/np.pi, eta_yumi, lw=3, label="Asymmetric Yumi")
ax.set_xlabel("Brace angle $\\theta_0$ (deg)", fontsize=16)
ax.set_ylabel("Static efficiency $\\eta_s$", fontsize=16)
ax.legend()
ax.set_title("(c)", loc="left", fontweight="bold")

# =====================================
# (d) ARROW VELOCITY VS MASS
# =====================================

masses = np.linspace(0.015, 0.05, 100)

E_self_tot = stored_energy(F_self_x, x_self)
E_rec_tot  = stored_energy(F_rec_x,  x_rec)
E_yumi_tot = stored_energy(F_yumi_x, x_yumi)

v_self = np.sqrt(2 * eta_d * E_self_tot / masses)
v_rec  = np.sqrt(2 * eta_d * E_rec_tot  / masses)
v_yumi = np.sqrt(2 * eta_d * E_yumi_tot / masses)

v_ref = np.sqrt(2 * eta_d * E_self_tot / m_arrow)

ax = axes[1, 1]
ax.plot(masses*1000, v_self / v_ref, lw=3, label="Self-Bow")
ax.plot(masses*1000, v_rec  / v_ref, lw=3, label="Composite Recurve")
ax.plot(masses*1000, v_yumi / v_ref, lw=3, label="Asymmetric Yumi")
ax.set_xlabel("Arrow mass (g)", fontsize=16)
ax.set_ylabel("Normalized launch velocity", fontsize=16)
ax.legend()
ax.set_title("(d)", loc="left", fontweight="bold")

# =====================================
# (e) MONTE CARLO
# =====================================

v_MC = monte_carlo_recurve_velocity(theta_max_rec)

ax = axes[2, 0]
ax.hist(v_MC, bins=30, density=True)
ax.set_xlabel("Arrow velocity $v$ (m/s)", fontsize=16)
ax.set_ylabel("Probability density", fontsize=16)
ax.set_title("(e)", loc="left", fontweight="bold")

# =====================================
# (f) MOMENT ARM
# =====================================

rc_self = ell_self * np.cos(ph_self - theta_self)
rc_rec  = ell_rec  * np.cos(ph_rec  - theta_rec) + \
          ell_r * np.exp(-alpha * np.maximum(theta_rec - theta0_rec, 0))
rc_yumi = 0.5 * (rc1 + rc2)

ax = axes[2, 1]
ax.plot(theta_self*180/np.pi, rc_self, lw=3, label="Self-Bow")
ax.plot(theta_rec*180/np.pi,  rc_rec,  lw=3, label="Composite Recurve")
ax.plot(theta1_yumi*180/np.pi, rc_yumi, lw=3, label="Asymmetric Yumi")
ax.set_xlabel("Limb rotation angle $\\theta$ (deg)", fontsize=16)
ax.set_ylabel("Effective moment arm $r_c$ (m)", fontsize=16)
ax.legend()
ax.set_title("(f)", loc="left", fontweight="bold")

fig.savefig("figures/Figure-composite.png", dpi=600)
plt.show()

# =====================================
# DYNAMIC EFFICIENCY SENSITIVITY
# =====================================

eta_vals = np.linspace(0.70, 0.90, 9)

v_eta = {
    "Self": np.sqrt(2 * eta_vals * E_self_tot / m_arrow),
    "Recurve": np.sqrt(2 * eta_vals * E_rec_tot / m_arrow),
    "Yumi": np.sqrt(2 * eta_vals * E_yumi_tot / m_arrow)
}

plt.figure(figsize=(4, 3))
plt.plot(eta_vals, v_eta["Self"], label="Self-Bow")
plt.plot(eta_vals, v_eta["Recurve"], label="Composite Recurve")
plt.plot(eta_vals, v_eta["Yumi"], label="Asymmetric Yumi")
plt.xlabel(r"Dynamic efficiency $\eta_d$")
plt.ylabel("Launch velocity (m/s)")
plt.legend()
plt.tight_layout()
plt.savefig("figures/sensitivity.png", dpi=600)
plt.show()
