# Standard imports
import numpy as np

class SemiPeriodicFourierSpline:

    def __init__(self, r: np.ndarray, angle: np.ndarray, F: np.ndarray, max_modenumber: int = 5):
        self._check_input(r, angle, F, max_modenumber)
        self.r = r
        self.modenumbers = np.arange(0, max_modenumber)
        self._compute_and_set_fourier_coefs(r, angle, F)
    
    def _check_input(self, r, angle, F, max_modenumber):
        if r.ndim != 1:
            raise ValueError("r must be a 1D array")
        if not is_monotonic(r):
            raise ValueError("r must be monotonically increasing/decreasing") 
        if angle.ndim != 1:
            raise ValueError("angle must be a 1D array")
        if not is_monotonic(angle):
            raise ValueError("angle must be monotonically increasing/decreasing")
        if F.ndim != 2:
            raise ValueError("F must be a 2D array")
        if F.shape != (len(r), len(angle)):
            raise ValueError("F must have shape (len(r), len(angle))")
        if not np.isreal(F).all():
            raise ValueError("F must be real")
        if max_modenumber < 0:
            raise ValueError("max_modenumber must be non-negative")
        pass

    def _compute_and_set_fourier_coefs(self, r, angle, F):
        from libneo import fourier_coefs_half
        self.fourier_coefs = np.zeros((len(r), len(self.modenumbers)), dtype=complex)
        for i, r in enumerate(r):
            self.fourier_coefs[i, :] = fourier_coefs_half(F[i, :-1], angle, self.modenumbers)

    def __call__(self, r, angle):
        coefs = self._get_fourier_coefs(r)
        return np.real(coefs.dot(np.exp(1.0j * np.outer(self.modenumbers, angle))))

    def _get_fourier_coefs(self, r):
        from scipy.interpolate import interp1d
        return interp1d(self.r, self.fourier_coefs, axis=0)(r)

def is_monotonic(x):
    return np.all(np.diff(x) > 0) or np.all(np.diff(x) < 0)