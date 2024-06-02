# Standard imports
import numpy as np

class SemiPeriodicFourierSpline:

    def __init__(self, r: np.ndarray, angle: np.ndarray, F: np.ndarray, max_modenumber: int=None):
        if angle.ndim == 1:
            angle = np.tile(angle, (len(r), 1))
        if max_modenumber is None:
            max_modenumber = (angle.shape[1]-1) // 2
        elif max_modenumber > (angle.shape[1]-1) // 2:
            raise SemiPeriodicFourierSpline_ExceedingNyquistLimit()
        self._check_input(r, angle, F, max_modenumber)
        self._r = r
        self.modenumbers = np.arange(0, max_modenumber+1)
        self._compute_and_set_fourier_coefs(r, angle, F)
    
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

    def _compute_and_set_fourier_coefs(self, r, angle, F):
        from libneo import get_half_fft
        self._fourier_coefs = np.zeros((len(r), len(self.modenumbers)), dtype=complex)
        for i, _ in enumerate(r):
            fm, _ = get_half_fft(F[i, :], angle[i,:])
            self._fourier_coefs[i, :] = fm[:len(self.modenumbers)]

    def __call__(self, r, angle, grid=True):
        coefs = self._get_fourier_coefs(r)
        if grid:
            return np.real(coefs.dot(np.exp(1.0j * np.outer(self.modenumbers, angle))))
        else:
            if isinstance(r, np.ndarray) and isinstance(angle, np.ndarray):
                if angle.ndim != 1:
                    raise ValueError("angle must be a 1D array")
                if r.shape != angle.shape:
                    raise ValueError("r and angle arrays must have the same shape if grid=False")
            return np.real(np.sum(coefs * (np.exp(1.0j * np.outer(angle, self.modenumbers))), axis=1))

    def _get_fourier_coefs(self, r):
        from scipy.interpolate import interp1d
        return interp1d(self._r, self._fourier_coefs, axis=0)(r)

def is_monotonic(x):
    return np.all(np.diff(x) > 0) or np.all(np.diff(x) < 0)

class SemiPeriodicFourierSpline_ExceedingNyquistLimit(Exception):
    def __init__(self):
        super().__init__()
    def __str__(self):
        return "The max modenumber must be less than or equal to (len(angle)-1) // 2"