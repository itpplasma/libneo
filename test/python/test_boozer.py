import pytest

import numpy as np
from libneo.boozer import (get_boozer_harmonics_divide_f_by_B0, get_boozer_harmonics,
    get_boozer_transform)
import matplotlib.pyplot as plt

def test_get_boozer_transform():
    stor = np.array([0.9, 0.6])
    nth = np.array([32])
    dth_of_thb, G_of_thb = get_boozer_transform(stor, nth)
    print(dth_of_thb)
    print(G_of_thb)
    

def test_get_boozer_harmonics_1D():
    nth = np.array([16])

    f = lambda spol, theta, phi: 2 *np.sin(theta) + 3 * np.cos(theta)

    stor = np.array([0.5, 0.6])
    dth_of_thb = []
    G_of_thb = []

    m0b = 16
    m = np.arange(-m0b, m0b+1)
    nph = 1
    plt.figure()
    for kth in nth:

        theta = np.linspace(0, 2*np.pi, kth, endpoint=False)
        phi = np.linspace(0, 2*np.pi, nph, endpoint=False)
        XX,YY = np.meshgrid(theta, phi)
        FF = f(stor[0], XX, YY)#np.ones_like(XX)
        f_fft = np.fft.fft2(FF)
        fmn_fft = f_fft / (kth * nph)

        print(kth)
        dth_of_thb.append(lambda thb: 0.0)
        G_of_thb.append(lambda thb: 0.0)

        res = get_boozer_harmonics(f, stor=stor, num_theta=kth, num_phi=nph, m0b=m0b,
            n=2, dth_of_thb=dth_of_thb, G_of_thb=G_of_thb) 

        plt.plot(m, np.abs(res[0]), label=str(kth))
        plt.plot(np.arange(kth), np.abs(fmn_fft[0,:]), label='fft, '+str(kth), ls = '--')
    plt.legend()
    plt.xlabel('m')
    plt.ylabel(r'abs($C_n$)')
    plt.show()


def test_get_boozer_harmonics_2D():
    nth = np.array([32])

    f = lambda spol, theta, phi: np.sin(5*theta) * np.sin(2 * phi) + np.sin(2 * theta) * np.sin(phi)

    stor = np.array([0.5, 0.6])
    dth_of_thb = []
    G_of_thb = []

    m0b = 24
    m = np.arange(-m0b, m0b + 1)
    nph = nth[0]
    plt.figure()
    for kth in nth:
        dth_of_thb.append(lambda thb: 0.0)
        G_of_thb.append(lambda thb: 0.5 * np.sin(thb))

        theta_b = np.linspace(0, 2 * np.pi, kth, endpoint=False)
        theta = theta_b + dth_of_thb[0](theta_b)
        phi_b = np.linspace(0, 2 * np.pi, nph, endpoint=False)
        phi = phi_b - G_of_thb[0](theta_b)

        XX, YY = np.meshgrid(theta, phi)
        FF = f(stor[0], XX, YY)
        f_fft = np.fft.fft2(FF)
        fmn_fft = f_fft / (kth * nph)

        print(kth)

        res = get_boozer_harmonics(
            f,
            stor=stor,
            num_theta=kth,
            num_phi=nph,
            m0b=m0b,
            n=2,
            dth_of_thb=dth_of_thb,
            G_of_thb=G_of_thb,
        )

        print(fmn_fft.shape)

        plt.plot(m, np.abs(res[0]), label=str(kth), lw=3)
        for n_mode in [2]: # check the toroidal modes corresponding to the ones in the defined function
            plt.plot(np.arange(kth), np.abs(fmn_fft[n_mode, :]), label="fft, n=" + str(n_mode), ls="--")
    
    plt.legend()
    plt.xlabel("m")
    plt.ylabel(r"abs($C_n$)")
    plt.show()


if __name__ == "__main__":
    pytest.main([__file__])
