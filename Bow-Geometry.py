import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve  

# =====================================
# BOW-SPECIFIC PARAMETERS (from Table 1)
# =====================================

# COMMON PARAMETERS
# Different max draws per bow type as per Table 1
X_MAX_SELF = 0.70    # meters
X_MAX_REC = 0.70     # meters  
X_MAX_YUMI = 0.90    # meters 
COS_EPS = 0.05


m_arrow = 0.025        # kg (25g arrow)
eta_d = 0.80           # dynamic efficiency

# ---------- SELF-BOW (Type A) ----------
ell_self = 0.75        # limb length (m)
y0_self = 0.10         # grip offset (m)
theta0_self = 20 * np.pi/180
k_self = 620.0         # rotational stiffness (Nm/rad)
L0_self = y0_self + ell_self * np.cos(theta0_self)
# ---------- RECURVE (Type B) ----------
ell_rec = 0.65         # limb length (m)
y0_rec = 0.10          # grip offset (m)
theta0_rec = 25 * np.pi/180  # brace + recurve offset
k_rec = 620.0          # rotational stiffness (Nm/rad)
L0_rec = y0_rec + ell_rec * np.cos(theta0_rec) 

# Recurve geometric parameters (from table)
ell_r = 0.03           # recurve length (m)
alpha = 2.0            # decay parameter 

# ---------- YUMI (Type C) ----------
ell1_yumi = 0.40       # shorter limb (m)
ell2_yumi = 0.88       # longer limb (m)
lambda_asym = 2.2      # limb length ratio
y0_yumi = 0.15         # grip offset (m)
theta0_yumi = 20 * np.pi/180
k1_yumi = 650.0        # Nm/rad (shorter limb)
k2_yumi = 590.0        # Nm/rad (longer limb) - average = 620
L0_yumi = y0_yumi + (ell1_yumi + ell2_yumi)/2 * np.cos(theta0_yumi)
ell_avg_yumi = (ell1_yumi + ell2_yumi) / 2
k_avg_yumi = (k1_yumi + k2_yumi) / 2



# =====================================
# UPDATED GEOMETRY FUNCTIONS (bow-specific)
# =====================================
def phi(theta, ell, y0, L0):
    """String angle for given bow geometry"""
    arg = (y0 + ell * np.cos(theta)) / L0
    arg = np.clip(arg, -1, 1)
    return np.arcsin(arg)

def draw_length(theta, ell, y0, L0, theta0):
    """Draw length for given bow geometry"""
    ph = phi(theta, ell, y0, L0)
    ph0 = phi(theta0, ell, y0, L0)
    return (ell * np.sin(theta) + L0 * np.cos(ph)) \
         - (ell * np.sin(theta0) + L0 * np.cos(ph0))


# =====================================
# FORCE LAWS (bow-specific)
# =====================================
def F_self(theta):
    ph = phi(theta, ell_self, y0_self, L0_self)
    return (2 * k_self * (theta - theta0_self) * np.cos(ph)) / \
           (ell_self * np.cos(ph - theta))

def F_recurve(theta):
    ph = phi(theta, ell_rec, y0_rec, L0_rec)
    rc = ell_rec * np.cos(ph - theta) + ell_r * np.exp(-alpha * np.maximum(theta - theta0_rec, 0))
    rc = np.maximum(rc, 1e-4)
    return (2 * k_rec * (theta - theta0_rec) * np.cos(ph)) / rc


def find_theta_max(ell, y0, L0, theta0, target_x):
    def f(theta):
        return draw_length(theta, ell, y0, L0, theta0) - target_x

    sol = fsolve(f, theta0 + 0.5, full_output=True)
    if sol[2] != 1:  # Check if solution converged
        # Try different initial guess
        sol = fsolve(f, theta0 + 1.0, full_output=True)
    return sol[0][0]






# Calculate max angles for X_MAX draw
theta_max_self = find_theta_max(ell_self, y0_self, L0_self, theta0_self, X_MAX_SELF)
theta_max_rec = find_theta_max(ell_rec, y0_rec, L0_rec, theta0_rec, X_MAX_REC)
theta_max_yumi = find_theta_max(ell_avg_yumi, y0_yumi, L0_yumi, theta0_yumi, X_MAX_YUMI)

# =====================================
# GENERATE ANGLE GRIDS FOR EACH BOW
# =====================================
theta_self = np.linspace(theta0_self + 1e-4, theta_max_self, 600)
theta_rec = np.linspace(theta0_rec + 1e-4, theta_max_rec, 600)
theta1_yumi = np.linspace(theta0_yumi + 1e-4, theta_max_yumi, 600)
theta2_yumi = np.zeros_like(theta1_yumi)

def yumi_tension_balance(theta2, theta1):
    ph1 = phi(theta1, ell1_yumi, y0_yumi, L0_yumi)
    ph2 = phi(theta2, ell2_yumi, y0_yumi, L0_yumi)

    rc1 = np.maximum(ell1_yumi * np.cos(ph1 - theta1), 1e-4)
    rc2 = np.maximum(ell2_yumi * np.cos(ph2 - theta2), 1e-4)

    return k1_yumi*(theta1 - theta0_yumi)/rc1 \
         - k2_yumi*(theta2 - theta0_yumi)/rc2

for i, th1 in enumerate(theta1_yumi):
    guess = theta2_yumi[i-1] if i > 0 else theta0_yumi + 0.1  # Start near brace angle
    theta2_yumi[i] = fsolve(
        yumi_tension_balance,
        guess,
        args=(th1,),
        xtol=1e-8,
        maxfev=1000
    )[0]



# Calculate draw lengths
x_self = draw_length(theta_self, ell_self, y0_self, L0_self,theta0_self)
x_rec = draw_length(theta_rec, ell_rec, y0_rec, L0_rec,theta0_rec)
x_yumi = (
    draw_length(theta1_yumi, ell1_yumi, y0_yumi, L0_yumi, theta0_yumi)
  + draw_length(theta2_yumi, ell2_yumi, y0_yumi, L0_yumi, theta0_yumi)
)


# Calculate string angles
ph_self = phi(theta_self, ell_self, y0_self, L0_self)
ph_rec  = phi(theta_rec,  ell_rec,  y0_rec,  L0_rec)
ph1 = phi(theta1_yumi, ell1_yumi, y0_yumi, L0_yumi)
ph2 = phi(theta2_yumi, ell2_yumi, y0_yumi, L0_yumi)


dx_dtheta_self = np.gradient(x_self, theta_self)
dx_dtheta_rec  = np.gradient(x_rec,  theta_rec)
dx_dtheta_yumi = np.gradient(x_yumi, theta1_yumi)

# Geometric admissibility masks
mask_self = (np.abs(np.cos(ph_self - theta_self)) > COS_EPS) & (dx_dtheta_self > 0)
mask_rec  = (np.abs(np.cos(ph_rec  - theta_rec )) > COS_EPS) & (dx_dtheta_rec  > 0)
mask_yumi = (
    (np.abs(np.cos(ph1 - theta1_yumi)) > COS_EPS)
    & (dx_dtheta_yumi > 0)
)

# Apply masks
theta_self, x_self = theta_self[mask_self], x_self[mask_self]
theta_rec,  x_rec  = theta_rec[mask_rec],  x_rec[mask_rec]
theta1_yumi = theta1_yumi[mask_yumi]
theta2_yumi = theta2_yumi[mask_yumi]
x_yumi      = x_yumi[mask_yumi]

# Calculate forces
F_self_x = F_self(theta_self)
F_rec_x = F_recurve(theta_rec)

# RECALCULATE string angles for filtered Yumi arrays
ph1 = phi(theta1_yumi, ell1_yumi, y0_yumi, L0_yumi)
ph2 = phi(theta2_yumi, ell2_yumi, y0_yumi, L0_yumi)

rc1 = np.maximum(ell1_yumi * np.cos(ph1 - theta1_yumi), 1e-4)
rc2 = np.maximum(ell2_yumi * np.cos(ph2 - theta2_yumi), 1e-4)

# Tension from balance (should be equal for both limbs)
T_yumi = k1_yumi * (theta1_yumi - theta0_yumi) / rc1
F_yumi_x = T_yumi * (np.cos(ph1) + np.cos(ph2))


dx = np.gradient(x_yumi, theta1_yumi)
mask = dx > 0

# =====================================
# Energy and efficiency utilities
# =====================================
def stored_energy(F, x):
    return np.trapezoid(F, x)

def static_efficiency(F, x):
    return stored_energy(F, x) / (np.max(F) * (x[-1] - x[0]))


def arrow_velocity(F, x):
    return np.sqrt(2 * eta_d * stored_energy(F, x) / m_arrow)





#==========Plots==============

fig, axes = plt.subplots(3, 2, figsize=(8.27, 11.69),constrained_layout=True)

# =====================================
# PANEL A: Force–draw curves
# =====================================
ax = axes[0,0]
ax.plot(x_self, F_self_x, lw=3, label="Self-Bow")
ax.plot(x_rec, F_rec_x, lw=3, label="Composite Recurve")
ax.plot(x_yumi, F_yumi_x, lw=3, label="Asymmetric Yumi")
ax.set_xlabel(r"Draw length $x$ (m)", fontsize=16)
ax.set_ylabel(r"Draw force $F$ (N)", fontsize=16)
ax.legend()
ax.set_title("(a)", loc="left", fontweight="bold")
# =====================================
# PANEL B: Stored energy vs draw
# =====================================
E_self = np.array([stored_energy(F_self_x[:i], x_self[:i]) for i in range(10, len(x_self))])
E_rec  = np.array([stored_energy(F_rec_x[:i],  x_rec[:i]) for i in range(10, len(x_rec))])
E_yumi = np.array([stored_energy(F_yumi_x[:i], x_yumi[:i]) for i in range(10, len(x_yumi))])


ax = axes[0,1]
ax.plot(x_self[10:], E_self, lw=3, label="Self-Bow")
ax.plot(x_rec[10:], E_rec,  lw=3, label="Composite Recurve")
ax.plot(x_yumi[10:], E_yumi, lw=3, label="Asymmetric Yumi")
ax.set_xlabel("Draw length $x$ (m)", fontsize=16)
ax.set_ylabel("Stored energy $U$ (J)", fontsize=16)
ax.legend()
ax.set_title("(b)", loc="left", fontweight="bold")


# =====================================
# PANEL C (REPLACED): Static efficiency vs brace angle
# =====================================
theta0_vals = np.linspace(10*np.pi/180, 40*np.pi/180, 50)

eta_self, eta_rec, eta_yumi = [], [], []

for th0 in theta0_vals:
    # Self-bow
    L0_temp = y0_self + ell_self * np.cos(th0) 
    def max_draw_eq(theta):
        return draw_length(theta, ell_self, y0_self, L0_temp, th0) - X_MAX_SELF


    sol = fsolve(max_draw_eq, th0 + 0.5, full_output=True, xtol=1e-8)
    if sol[2] != 1:  # If didn't converge
        # Try different initial guesses
        for guess in [th0 + 0.3, th0 + 0.7, th0 + 1.0]:
            sol = fsolve(max_draw_eq, guess, full_output=True, xtol=1e-8)
            if sol[2] == 1:
                break
    theta_max_temp = sol[0][0]
    
    theta_temp = np.linspace(th0 + 1e-4, theta_max_temp, 400)

    ph_temp = phi(theta_temp, ell_self, y0_self, L0_temp)
    x_temp = draw_length(theta_temp, ell_self, y0_self, L0_temp, th0)
    F_s = (2 * k_self * (theta_temp - th0) * np.cos(ph_temp)) / \
          (ell_self * np.cos(ph_temp - theta_temp))
    eta_self.append(static_efficiency(F_s, x_temp))

    # Recurve
    L0_temp = y0_rec + ell_rec * np.cos(th0) 
    theta_max_temp = find_theta_max(ell_rec, y0_rec, L0_temp, th0, X_MAX_REC)
    theta_temp = np.linspace(th0 + 1e-4, theta_max_temp, 400)
    ph_temp = phi(theta_temp, ell_rec, y0_rec, L0_temp)
    x_temp = draw_length(theta_temp, ell_rec, y0_rec, L0_temp, th0)

    #x_temp = ell_rec * np.sin(theta_temp) + L0_temp * np.cos(ph_temp)
    F_r = (2 * k_rec * (theta_temp - th0) * np.cos(ph_temp)) / \
          (ell_rec * np.cos(ph_temp - theta_temp) +
           ell_r * np.exp(-alpha * (theta_temp - th0)))
    eta_rec.append(static_efficiency(F_r, x_temp))
    
    # Yumi 
    ell_avg = (ell1_yumi + ell2_yumi) / 2
    k_avg = (k1_yumi + k2_yumi) / 2
    L0_temp = y0_yumi + ell_avg * np.cos(th0)
    theta_max_temp = find_theta_max(ell_avg, y0_yumi, L0_temp, th0, X_MAX_YUMI)
    theta_temp = np.linspace(th0 + 1e-4, theta_max_temp, 400)
    ph_temp = phi(theta_temp, ell_avg, y0_yumi, L0_temp)
    x_temp = draw_length(theta_temp, ell_avg, y0_yumi, L0_temp, th0)
    
    # CORRECTION: For asymmetric, F = 2T cos φ only if limbs are symmetric
    # With average parameters, we use the symmetric formula
    F_y = (2 * k_avg * (theta_temp - th0) * np.cos(ph_temp)) / \
          (ell_avg * np.cos(ph_temp - theta_temp))
    eta_yumi.append(static_efficiency(F_y, x_temp))
ax = axes[1,0]
ax.plot(theta0_vals*180/np.pi, eta_self, lw=3, label="Self-Bow")
ax.plot(theta0_vals*180/np.pi, eta_rec, lw=3, label="Composite Recurve")
ax.plot(theta0_vals*180/np.pi, eta_yumi, lw=3, label="Asymmetric Yumi")
ax.set_xlabel("Brace angle $\\theta_0$ (deg)", fontsize=16)
ax.set_ylabel("Static efficiency $\\eta_s$", fontsize=16)
ax.legend()
ax.set_title("(c)", loc="left", fontweight="bold")



# =====================================
# PANEL D: Arrow velocity vs arrow mass
# =====================================
masses = np.linspace(0.015, 0.05, 100)
E_self_tot = stored_energy(F_self_x, x_self)
E_rec_tot  = stored_energy(F_rec_x,  x_rec)
E_yumi_tot = stored_energy(F_yumi_x, x_yumi)

v_self = np.sqrt(2 * eta_d * E_self_tot / masses)
v_rec  = np.sqrt(2 * eta_d * E_rec_tot  / masses)
v_yumi = np.sqrt(2 * eta_d * E_yumi_tot / masses)

v_ref = np.sqrt(2 * eta_d * E_self_tot / m_arrow)


ax = axes[1,1]
ax.plot(masses*1000, v_self / v_ref, lw=3, label="Self-Bow")
ax.plot(masses*1000, v_rec / v_ref,  lw=3, label="Composite Recurve")
ax.plot(masses*1000, v_yumi/ v_ref, lw=3, label="Asymmetric Yumi")
ax.set_xlabel("Arrow mass (g)", fontsize=16)
ax.set_ylabel("Normalized launch velocity", fontsize=16)
ax.legend()
ax.set_title("(d)", loc="left", fontweight="bold")


# =====================================
# PANEL E: Monte Carlo sensitivity
# =====================================
np.random.seed(0)
N = 1000
k_var = np.random.normal(1.0, 0.15, N)
ell_var = np.random.normal(1.0, 0.05, N)

v_MC = []

# Vary k and recalculate force curve
for kf, lf in zip(k_var, ell_var):
    # Create temporary parameters
    k_temp = kf * k_rec
    ell_temp = lf * ell_rec
    
    # Recalculate force with new parameters
    def F_MC_func(theta):
        ph = phi(theta, ell_temp, y0_rec, L0_rec)
        rc = ell_temp * np.cos(ph - theta) + ell_r * np.exp(-alpha * (theta - theta0_rec))
        return (2 * k_temp * (theta - theta0_rec) * np.cos(ph)) / rc
    
    # Need to recalculate x too
    theta_MC = np.linspace(theta0_rec + 1e-4, theta_max_rec, 600)
    x_MC = draw_length(theta_MC, ell_temp, y0_rec, L0_rec, theta0_rec)
    F_MC = F_MC_func(theta_MC)

    dx = np.gradient(x_MC, theta_MC)
    mask = dx > 0
    x_MC = x_MC[mask]
    F_MC = F_MC[mask]

    v_MC.append(arrow_velocity(F_MC, x_MC))

ax = axes[2,0]
ax.hist(v_MC, bins=30, density=True)
ax.set_xlabel("Arrow velocity $v$ (m/s)", fontsize=16)
ax.set_ylabel("Probability density", fontsize=16)
ax.set_title("(e)", loc="left", fontweight="bold")


# =====================================
# PANEL F: Effective moment arm vs draw angle
# =====================================
# Self-bow moment arm
ph_self = phi(theta_self, ell_self, y0_self, L0_self)
rc_self = ell_self * np.cos(ph_self - theta_self)

# Recurve moment arm
ph_rec = phi(theta_rec, ell_rec, y0_rec, L0_rec)
rc_rec = ell_rec * np.cos(ph_rec - theta_rec) + \
         ell_r * np.exp(-alpha * np.maximum(theta_rec - theta0_rec, 0))


# Yumi moment arm (average)
rc_yumi = 0.5 * (rc1 + rc2)



ax = axes[2,1]
ax.plot(theta_self*180/np.pi, rc_self, lw=3, label="Self-Bow")
ax.plot(theta_rec*180/np.pi, rc_rec, lw=3, label="Composite Recurve")
ax.plot(theta1_yumi*180/np.pi, rc_yumi, lw=3, label="Asymmetric Yumi")
ax.set_xlabel("Limb rotation angle $\\theta$ (deg)", fontsize=16)
ax.set_ylabel("Effective moment arm $r_c$ (m)", fontsize=16)
ax.legend()
ax.set_title("(f)", loc="left", fontweight="bold")


fig.set_constrained_layout_pads(wspace=0.15, hspace=0.15)
fig.savefig("Figure-composite.png", dpi=600)
plt.show()

#========== Sensitivity plots==============
# =====================================
# Dynamic efficiency sensitivity test
# =====================================
eta_vals = np.linspace(0.70, 0.90, 9)

v_eta = {
    "Self": [],
    "Recurve": [],
    "Yumi": []
}

for eta in eta_vals:
    v_eta["Self"].append(np.sqrt(2 * eta * E_self_tot / m_arrow))
    v_eta["Recurve"].append(np.sqrt(2 * eta * E_rec_tot  / m_arrow))
    v_eta["Yumi"].append(np.sqrt(2 * eta * E_yumi_tot / m_arrow))

v_eta = {k: np.array(v) for k, v in v_eta.items()}

plt.figure(figsize=(4,3))
plt.plot(eta_vals, v_eta["Self"], label="Self-Bow")
plt.plot(eta_vals, v_eta["Recurve"], label="Composite Recurve")
plt.plot(eta_vals, v_eta["Yumi"], label="Asymmetric Yumi")
plt.xlabel(r"Dynamic efficiency $\eta_d$")
plt.ylabel("Launch velocity (m/s)")
plt.legend()
plt.tight_layout()
plt.savefig('sensitivity.png', dpi=600)
plt.show()
