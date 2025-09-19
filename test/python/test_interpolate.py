import numpy as np
import pytest
from libneo.interpolate import Interpolate

# XFAIL until a proper Fortran-backed Python binding for Interpolate is shipped
pytestmark = pytest.mark.xfail(reason="Python binding for Fortran 'interpolate' not built in CI; see issue #127", strict=False)

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
