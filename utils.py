import numpy as np
from scipy.optimize import fsolve  

# =====================================
# ENERGY AND PERFORMANCE UTILITIES
# =====================================

def stored_energy(F, x):
    """
    Compute stored elastic energy using numerical integration.

    Parameters
    ----------
    F : array_like
        Draw force array (N)
    x : array_like
        Corresponding draw length array (m)

    Returns
    -------
    float
        Stored energy (J)
    """
    return np.trapz(F, x)


def static_efficiency(F, x):
    """
    Compute static efficiency defined as:
    eta_s = U / [F_max * (x_max - x_min)]

    Parameters
    ----------
    F : array_like
        Draw force array (N)
    x : array_like
        Draw length array (m)

    Returns
    -------
    float
        Static efficiency (dimensionless)
    """
    return stored_energy(F, x) / (np.max(F) * (x[-1] - x[0]))


def arrow_velocity(F, x, eta_d, m_arrow):
    """
    Compute arrow launch velocity assuming fixed dynamic efficiency.

    Parameters
    ----------
    F : array_like
        Draw force array (N)
    x : array_like
        Draw length array (m)
    eta_d : float
        Dynamic efficiency
    m_arrow : float
        Arrow mass (kg)

    Returns
    -------
    float
        Arrow launch velocity (m/s)
    """
    return np.sqrt(2 * eta_d * stored_energy(F, x) / m_arrow)
