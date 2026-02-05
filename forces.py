import numpy as np
from geometry import phi
from parameters import *


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


