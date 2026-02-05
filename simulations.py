import numpy as np
from scipy.optimize import fsolve  
from geometry import phi, draw_length
from parameters import *
import warnings

# =====================================
# GENERIC SOLVER UTILITIES
# =====================================

def find_theta_max(ell, y0, L0, theta0, target_x):
    """
    Solve for maximum limb rotation angle corresponding
    to a target draw length.
    """
    def f(theta):
        return draw_length(theta, ell, y0, L0, theta0) - target_x

    sol = fsolve(f, theta0 + 0.5, full_output=True)
    if sol[2] != 1:
        sol = fsolve(f, theta0 + 1.0, full_output=True)

    return sol[0][0]

# =====================================
# YUMI ASYMMETRIC LIMB COUPLING
# =====================================

def yumi_tension_balance(theta2, theta1):
    """
    Enforces equal string tension in asymmetric Yumi limbs.
    """
    ph1 = phi(theta1, ell1_yumi, y0_yumi, L0_yumi)
    ph2 = phi(theta2, ell2_yumi, y0_yumi, L0_yumi)

    rc1 = np.maximum(ell1_yumi * np.cos(ph1 - theta1), 1e-4)
    rc2 = np.maximum(ell2_yumi * np.cos(ph2 - theta2), 1e-4)

    return (
        k1_yumi * (theta1 - theta0_yumi) / rc1
        - k2_yumi * (theta2 - theta0_yumi) / rc2
    )


def solve_yumi_secondary_angle(theta1_array):
    """
    Solve secondary limb angle theta2 for each theta1
    using force balance.
    """
    theta2 = np.zeros_like(theta1_array)

    for i, th1 in enumerate(theta1_array):
        guess = theta2[i - 1] if i > 0 else theta0_yumi + 0.1
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            theta2[i] = fsolve(
                yumi_tension_balance,
                guess,
                args=(th1,),
                xtol=1e-8,
                maxfev=2000
            )[0]


    return theta2




# =====================================
# MONTE CARLO SENSITIVITY (RECURVE)
# =====================================

def monte_carlo_recurve_velocity(
    theta_max,
    N=1000,
    seed=0
):
    """
    Monte Carlo sensitivity analysis for recurve bow
    with stiffness and limb-length perturbations.
    """
    np.random.seed(seed)

    k_var = np.random.normal(1.0, 0.15, N)
    ell_var = np.random.normal(1.0, 0.05, N)

    v_MC = []

    for kf, lf in zip(k_var, ell_var):
        k_temp = kf * k_rec
        ell_temp = lf * ell_rec

        def F_MC(theta):
            ph = phi(theta, ell_temp, y0_rec, L0_rec)
            rc = ell_temp * np.cos(ph - theta) + \
                 ell_r * np.exp(-alpha * (theta - theta0_rec))
            return (2 * k_temp * (theta - theta0_rec) * np.cos(ph)) / rc

        theta_MC = np.linspace(theta0_rec + 1e-4, theta_max, 600)
        x_MC = draw_length(theta_MC, ell_temp, y0_rec, L0_rec, theta0_rec)
        F_MC_vals = F_MC(theta_MC)

        dx = np.gradient(x_MC, theta_MC)
        mask = dx > 0

        x_MC = x_MC[mask]
        F_MC_vals = F_MC_vals[mask]

        from utils import arrow_velocity
        v_MC.append(
            arrow_velocity(F_MC_vals, x_MC, eta_d, m_arrow)
        )

    return np.array(v_MC)



