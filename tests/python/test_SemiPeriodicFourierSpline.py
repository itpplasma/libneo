# %% Standard imports
import numpy as np

# modules to test
from libneo import SemiPeriodicFourierSpline

r = np.linspace(0, 1, 5)
angle = np.linspace(0, 2*np.pi, 10)

def test_init():
    Angle, R = np.meshgrid(angle, r)
    F = R*np.cos(Angle)
    spline = SemiPeriodicFourierSpline(r, angle, F)

if __name__ == '__main__':
    test_init()