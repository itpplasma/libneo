import pytest

import numpy as np
from libneo.boozer import get_boozer_harmonics_divide_f_by_B0

def test_get_boozer_harmonics_divide_f_by_B0():
    nth = 8

    f = lambda spol, theta, phi: sin(theta + 2*phi)

    stor = np.array([0.5, 0.6])
    dth_of_thb = []
    G_of_thb = []

    for kth in range(nth):
        dth_of_thb.append(lambda thb: 0.0)
        G_of_thb.append(lambda thb: 0.0)

    get_boozer_harmonics_divide_f_by_B0(f, stor=stor, nth=nth, nph=16, m0b=24,
        n=2, dth_of_thb=dth_of_thb, G_of_thb=G_of_thb)

    print("YEYYY")

if __name__ == "__main__":
    pytest.main([__file__])
