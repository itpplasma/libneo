import numpy as np
from scipy.interpolate import CubicSpline, PPoly

def fourier_coefs_half(f, x, m):
    """
    x ranges over the full period, i.e. x[0] = x[-1] + period.
    length of f is len(x) - 1 and skips value for x[-1].
    reason: avoid double counting in Fourier integral
    while still allowing to get intervals by diff(x).
    """
    num_x = len(x)
    dx = np.diff(x)
    expfac = np.exp(-1.0j * np.outer(m, x[:-1]))
    fm = expfac.dot(f * dx) / (2 * np.pi)
    fm[1:] = 2.0 * fm[1:]  # Half-sided Fourier series
    return fm

def fourier_coefs_full(f, x, m):
    num_x = len(x)
    dx = np.diff(x)
    expfac = np.exp(-1.0j * np.outer(m, x[:-1]))
    fm = expfac.dot(f * dx) / (2 * np.pi)
    return fm
