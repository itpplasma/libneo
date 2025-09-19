import numpy as np
from libneo.interpolate import Interpolate

def test_construct_splines_1d():
    import libneo
    x_min = 0.0
    x_max = 1.0
    y = np.linspace(0, 1, 10)
    order = 5
    periodic = False
    spl = Interpolate.construct_splines_1d(x_min, x_max, y, order, periodic)
    assert type(spl) == Interpolate.SplineData1D

if __name__ == '__main__':
    test_construct_splines_1d()
