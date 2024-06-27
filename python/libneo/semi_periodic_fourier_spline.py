# Standard imports
import numpy as np
from scipy.interpolate import CubicSpline

class SemiPeriodicFourierSpline:

    def __init__(self, r: np.ndarray, angle: np.ndarray, F: np.ndarray, max_modenumber: int=None):
        if angle.ndim == 1:
            angle = np.tile(angle, (len(r), 1))
        if max_modenumber is None:
            max_modenumber = (angle.shape[1]-1) // 2
        elif max_modenumber > (angle.shape[1]-1) // 2:
            raise SemiPeriodicFourierSpline_ExceedingNyquistLimit()
        self._check_input(r, angle, F, max_modenumber)
        self._modenumbers = np.arange(0, max_modenumber+1)
        self._fourier_coefs = self._compute_fourier_coefs(r, angle, F)
        self._fourier_coefs_splines = CubicSpline(r, self._fourier_coefs, axis=0)
    
    def _check_input(self, r, angle, F, max_modenumber):
        if r.ndim != 1:
            raise ValueError("r must be a 1D array")
        if not is_monotonic(r):
            raise ValueError("r must be monotonically increasing/decreasing") 
        if angle.ndim > 1 and (angle.shape[0] != r.shape[0]):
            raise ValueError("angle must be a 1D array or have the same 0th dimension as r")
        if not is_monotonic(angle):
            raise ValueError("angle must be monotonically increasing/decreasing")
        if F.ndim != 2:
            raise ValueError("F must be a 2D array")
        if (F.shape != (len(r), len(angle))) and (F.shape != (len(r), angle.shape[1])):
            raise ValueError("F must have shape (len(r), len(angle)) or (len(r), angle.shape[1])")
        if not np.isreal(F).all():
            raise ValueError("F must be real")
        if max_modenumber < 0:
            raise ValueError("max_modenumber must be non-negative")
        pass

    def _compute_fourier_coefs(self, r, angle, F):
        from libneo import get_half_fft
        fourier_coefs = np.zeros((len(r), len(self._modenumbers)), dtype=complex)
        for i, _ in enumerate(r):
            fm, _ = get_half_fft(F[i, :], angle[i,:])
            fourier_coefs[i, :] = fm[:len(self._modenumbers)]
        return fourier_coefs

    def __call__(self, r, angle, grid=True):
        coefs = self._fourier_coefs_splines(r)
        if grid:
            return np.real(coefs.dot(np.exp(1.0j * np.outer(self._modenumbers, angle))))
        else:
            if isinstance(r, np.ndarray) and isinstance(angle, np.ndarray):
                if angle.ndim != 1:
                    raise ValueError("angle must be a 1D array")
                if r.shape != angle.shape:
                    raise ValueError("r and angle arrays must have the same shape if grid=False")
            return np.real(np.sum(coefs * (np.exp(1.0j * np.outer(angle, self._modenumbers))), axis=1))

def is_monotonic(x):
    return np.all(np.diff(x) > 0) or np.all(np.diff(x) < 0)

class SemiPeriodicFourierSpline_ExceedingNyquistLimit(Exception):
    def __init__(self):
        super().__init__()
    def __str__(self):
        return "The max modenumber must be less than or equal to (len(angle)-1) // 2"