import numpy as np


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


