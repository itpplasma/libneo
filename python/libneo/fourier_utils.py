# %%Standard imports
import numpy as np
from scipy.interpolate import CubicSpline

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

def get_half_fft(f, x):
    from scipy.interpolate import CubicSpline
    if np.iscomplexobj(f):
        raise ValueError("Function has to be real valued!")
    xi_equi = np.linspace(x[0], x[-1], len(x))[:-1] # skip last point that is repeated as
    sample_rate = (xi_equi[1] - xi_equi[0])         # DFT sum goes for xi only until i=N-1
    fi_equi = CubicSpline(x, f, extrapolate='periodic')(xi_equi)
    fm = np.fft.rfft(fi_equi)/len(fi_equi) # normalisation absorbed into coefs
    if len(xi_equi) % 2 == 0:
        fm[1:-1] *= 2.0 # for even samples, the 0th and highest mode dont have hermitian counterparts
    else:
        fm[1:] *= 2.0 # for odd samples, the only the 0th mode has no hermitian counterpart
    m = np.fft.rfftfreq(len(xi_equi), d=sample_rate)*(x[-1]-x[0])
    fm *= np.exp(-1.0j * m * x[0]) # shift to interval [0,2pi]
    return fm, m