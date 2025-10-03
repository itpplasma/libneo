from pathlib import Path
import os

import pytest
import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from _efit_to_boozer import efit_to_boozer
from efit_to_boozer.boozer import (get_boozer_harmonics_divide_f_by_B0, get_boozer_harmonics,
    get_boozer_transform, get_B0_of_s_theta_boozer, get_boozer_harmonics_divide_f_by_B0_1D,
    get_magnetic_axis)


DEBUG = os.environ.get("LIBNEO_DEBUG_BOOZER") == "1"


def _debug(msg: str) -> None:
    if DEBUG:
        print(f"[boozer-test] {msg}", flush=True)


@pytest.fixture(autouse=True)
def _set_test_cwd(monkeypatch):
    target = Path(__file__).parent
    monkeypatch.chdir(target)
    _debug(f"cwd set to {Path.cwd()}")
    repo_root = target.resolve().parents[1]
    convex_src = (target / '../../test/resources/convexwall.dat').resolve()
    _debug(f"convexwall source: {convex_src}")
    _debug(f"convexwall source exists: {convex_src.exists()}")
    _debug(f"test.geqdsk present: {(target / 'test.geqdsk').exists()}")
    _debug(f"field_divB0.inp present: {(target / 'field_divB0.inp').exists()}")


def _figure_path(name: str) -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    target = repo_root / "build" / "test"
    target.mkdir(parents=True, exist_ok=True)
    return target / name


def test_get_boozer_transform():
    _debug("start test_get_boozer_transform")
    stor = np.array([0.9, 0.6])
    nth = 32
    _debug("calling get_boozer_transform")
    (r_of_thb, dth_of_thb, G_of_thb, theta_boozers,
     theta_geoms, theta_symflux, B0) = get_boozer_transform(stor, nth)
    assert len(r_of_thb) == stor.size - 1
    assert len(dth_of_thb) == stor.size - 1
    assert len(G_of_thb) == stor.size - 1
    assert theta_boozers.shape == (stor.size - 1, 2 * nth + 1)
    assert theta_geoms.shape == (stor.size - 1, 2 * nth + 1)
    assert theta_symflux.shape == (2 * nth,)
    assert len(B0) == stor.size - 1
    _debug("done test_get_boozer_transform")
def test_get_boozer_harmonics_1D():
    _debug("start test_get_boozer_harmonics_1D")
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

        dth_of_thb.append(lambda thb: 0.0)
        G_of_thb.append(lambda thb: 0.0)

        _debug(f"computing harmonics for kth={kth}")
        res = get_boozer_harmonics(f, stor=stor, num_theta=kth, num_phi=nph, m0b=m0b,
            n=2, dth_of_thb=dth_of_thb, G_of_thb=G_of_thb)

        plt.plot(m, np.abs(res[0]), label=str(kth))
        plt.plot(np.arange(kth), np.abs(fmn_fft[0,:]), label='fft, '+str(kth), ls = '--')
    plt.legend()
    plt.xlabel('m')
    plt.ylabel(r'abs($C_n$)')
    plt.savefig(_figure_path('test_get_boozer_harmonics_1D.png'))
    _debug("done test_get_boozer_harmonics_1D")


def test_get_boozer_harmonics_2D():
    _debug("start test_get_boozer_harmonics_2D")
    nth = np.array([32])

    f = lambda spol, theta, phi: np.sin(5*theta) * np.sin(2 * phi) + 0.5*np.sin(2 * theta) * np.sin(2*phi)

    stor = np.array([0.5, 0.6])
    dth_of_thb = []
    G_of_thb = []

    m0b = 24
    m = np.arange(-m0b, m0b + 1)
    nph = nth[0]

    dth_of_thb.append(lambda thb: 0.1 * np.cos(thb))
    G_of_thb.append(lambda thb: 0.5 * np.sin(thb))

    plt.figure()

    for kth in nth:
        theta_b = np.linspace(0, 2 * np.pi, kth, endpoint=False)
        theta = theta_b + dth_of_thb[0](theta_b)
        phi_b = np.linspace(0, 2 * np.pi, nph, endpoint=False)
        phi = np.zeros((kth, nph))
        for i, phib in enumerate(phi_b):
            phi[:,i] = phib - G_of_thb[0](theta_b)

        THTH = np.tile(theta, (nph,1))
        PHIPHI = phi
        FF = f(stor[0], THTH, PHIPHI.T)
        f_fft = np.fft.fft2(FF)
        fmn_fft = f_fft / (kth * nph)

        _debug(f"computing 2D harmonics for kth={kth}")
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
    plt.savefig(_figure_path('test_get_boozer_harmonics_2D.png'))
    _debug("done test_get_boozer_harmonics_2D")


def test_get_boozer_harmonics_divide_f_by_B0_1D():
    _debug("start test_get_boozer_harmonics_divide_f_by_B0_1D")
    nth = np.array([64])

    n = 2
    f = lambda spol, theta, phi: (2 *np.sin(theta) + 3 * np.cos(3*theta)) * 10**4

    stor = np.array([0.5, 0.6])
    dth_of_thb = []
    G_of_thb = []

    m0b = 16
    m = np.arange(-m0b, m0b+1)
    nph = 1
    dth_of_thb.append(lambda thb: 0.3 * np.cos(thb))
    G_of_thb.append(lambda thb: 0.5 * np.sin(thb))
    stor_ind = 0

    plt.figure()

    for kth in nth:
        B0 = get_B0_of_s_theta_boozer(stor, kth)

        theta_b = np.linspace(0, 2 * np.pi, kth*2+1, endpoint=False)
        theta = theta_b + dth_of_thb[stor_ind](theta_b)

        B0_arr = np.array(B0[stor_ind](theta_b))

        phi_b = np.linspace(0, 2 * np.pi, nph, endpoint=False)
        phi = np.zeros((kth*2+1, nph))
        for i, phib in enumerate(phi_b):
            phi[:,i] = phib - G_of_thb[stor_ind](theta_b)

        THTH = theta
        PHIPHI = phi
        FF = f(stor[stor_ind], THTH, 0.0) / B0_arr * np.exp(-1j* n * G_of_thb[stor_ind](theta_b))
        f_fft = np.fft.fft(FF)
        fmn_fft = f_fft / ((2*kth+1) * nph)
        fmn_fft = np.fft.fftshift(fmn_fft)

        _debug(f"computing divide_f_by_B0_1D for kth={kth}")
        res = get_boozer_harmonics_divide_f_by_B0_1D(f, stor=stor, num_theta=kth, num_phi=nph, m0b=m0b,
            n=2, dth_of_thb=dth_of_thb, G_of_thb=G_of_thb)

        ind = int(((2*kth+1) - (2*m0b + 1))/ 2 )

        plt.plot(m, np.abs(res[0]), label=str(kth))
        plt.plot(m, np.abs(fmn_fft[ind:(2*kth+1)-ind]), label='fft, '+str(kth), ls = '--')
        #assert np.isclose(np.abs(res[0]), np.abs(fmn_fft[0,ind:(2*kth+1)-ind]), atol=1e-3).all()
    plt.legend()
    plt.xlabel('m')
    plt.ylabel(r'abs($C_n$)')
    plt.savefig(_figure_path('test_get_boozer_harmonics_divide_f_by_B0_1D.png'))
    _debug("done test_get_boozer_harmonics_divide_f_by_B0_1D")


def test_get_boozer_harmonics_divide_f_by_B0_2D():
    _debug("start test_get_boozer_harmonics_divide_f_by_B0_2D")
    nth = np.array([64])

    #f = lambda spol, theta, phi: 2*np.sin(theta) * np.sin(2 * phi) + 0.5*np.sin(3 * theta) * np.sin(2*phi)
    f = lambda spol, theta, phi: (2 *np.sin(theta) + 3 * np.cos(3*theta)) * 10**4 * np.sin(2 * phi)

    stor = np.array([0.5, 0.6])
    dth_of_thb = []
    G_of_thb = []

    m0b = 18
    m = np.arange(-m0b, m0b + 1)
    nph = nth[0]

    dth_of_thb.append(lambda thb: 0.3 * np.cos(thb))
    G_of_thb.append(lambda thb: 0.5 * np.sin(thb))

    plt.figure()

    stor_ind = 0

    for kth in nth:
        B0 = get_B0_of_s_theta_boozer(stor, kth)

        theta_b = np.linspace(0, 2 * np.pi, 2*kth+1, endpoint=False)
        theta = theta_b + dth_of_thb[stor_ind](theta_b)

        B0_arr = B0[stor_ind](theta_b)

        phi_b = np.linspace(0, 2 * np.pi, nph, endpoint=False)
        phi = np.zeros((2*kth+1, nph))
        for i, phib in enumerate(phi_b):
            phi[:,i] = phib - G_of_thb[stor_ind](theta_b)

        THTH = np.tile(theta, (nph,1))
        PHIPHI = phi
        FF = f(stor[0], THTH, PHIPHI.T) / B0_arr
        f_fft = np.fft.fft2(FF)
        fmn_fft = f_fft / ((kth*2+1) * nph)
        fmn_fft = np.fft.fftshift(fmn_fft, axes=1)

        _debug(f"computing divide_f_by_B0_2D for kth={kth}")
        res = get_boozer_harmonics_divide_f_by_B0(
            f,
            stor=stor,
            num_theta=kth,
            num_phi=nph,
            m0b=m0b,
            n=2,
            dth_of_thb=dth_of_thb,
            G_of_thb=G_of_thb,
        )

        ind = int(((2*kth+1) - (2*m0b + 1))/ 2 )

        plt.plot(m, np.abs(res[0]), label=str(kth), lw=3)
        for n_mode in [2]: # check the toroidal modes corresponding to the ones in the defined function
            plt.plot(m, np.abs(fmn_fft[n_mode, ind:(2*kth+1)-ind]), label="fft, n=" + str(n_mode), ls="--")

    plt.legend()
    plt.xlabel("m")
    plt.ylabel(r"abs($C_{mn}$)")
    plt.savefig(_figure_path('test_get_boozer_harmonics_divide_f_by_B0_2D.png'))
    _debug("done test_get_boozer_harmonics_divide_f_by_B0_2D")


@pytest.mark.skip(reason="Not implemented currently")
def test_get_boozer_harmonics_divide_f_by_B0_1D_fft():
    from efit_to_boozer import get_boozer_harmonics_divide_f_by_B0_1D_fft
    nth = np.array([64])

    n = 2
    f = lambda spol, theta: (2 *np.sin(theta) + 3 * np.cos(3*theta)) * 10**4

    stor = np.array([0.5, 0.6])
    dth_of_thb = []
    G_of_thb = []

    m0b = 16
    m = np.arange(-m0b, m0b+1)
    nph = 1
    dth_of_thb.append(lambda thb: 0.3 * np.cos(thb))
    G_of_thb.append(lambda thb: 0.5 * np.sin(thb))
    stor_ind = 0

    plt.figure()

    for kth in nth:
        B0 = get_B0_of_s_theta_boozer(stor, kth)

        theta_b = np.linspace(0, 2 * np.pi, kth, endpoint=False)
        theta = theta_b + dth_of_thb[stor_ind](theta_b)

        B0_arr = np.array(B0[stor_ind](theta_b))

        phi_b = np.linspace(0, 2 * np.pi, nph, endpoint=False)
        phi = np.zeros((kth, nph))
        for i, phib in enumerate(phi_b):
            phi[:,i] = phib - G_of_thb[stor_ind](theta_b)

        THTH = theta
        PHIPHI = phi
        FF = f(stor[stor_ind], THTH) / B0_arr * np.exp(-1j* n * G_of_thb[stor_ind](theta_b))
        f_fft = np.fft.fft(FF)
        fmn_fft = f_fft / ((kth) * nph)
        fmn_fft = np.fft.fftshift(fmn_fft)

        res = get_boozer_harmonics_divide_f_by_B0_1D_fft(f, stor=stor, num_theta=kth,
            n=2, dth_of_thb=dth_of_thb, G_of_thb=G_of_thb)
        res_modes = np.concatenate((np.arange(0, kth/2), np.arange(-kth/2, 0)))
        print(res_modes)

        ind = int(((kth) - (2*m0b + 1))/ 2 )

        plt.scatter(res_modes, np.abs(res[0,:]), label=str(kth))
        plt.plot(np.arange(-nth/2,nth/2), np.abs(fmn_fft), label='fft, '+str(kth), ls = '--')
        assert np.isclose(np.abs(res[0]), np.abs(fmn_fft[0,ind:(2*kth+1)-ind]), atol=1e-3).all()
    plt.legend()
    plt.xlabel('m')
    plt.ylabel(r'abs($C_n$)')
    plt.savefig('test_get_boozer_harmonics_divide_f_by_B0_1D_fft.png')


def test_get_B0_of_s_theta_boozer():
    _debug("start test_get_B0_of_s_theta_boozer")
    stor = np.array([0.3, 0.6])
    nth = 32
    theta = np.linspace(0, 2*np.pi, nth)
    _debug("calling get_B0_of_s_theta_boozer")
    B0 = get_B0_of_s_theta_boozer(stor, nth)
    plt.figure()
    plt.plot(theta, B0[0](theta))
    plt.ylabel(r"$B_0$ [G]")
    plt.xlabel(r"$\vartheta$ [rad]")
    plt.title(f"stor = {stor[0]}, nth = {nth}")
    plt.savefig(_figure_path('test_get_B0_of_s_theta_boozer.png'))
    _debug("done test_get_B0_of_s_theta_boozer")

def test_get_magnetic_axis():
    _debug("start test_get_magnetic_axis")
    efit_to_boozer.init()
    try:
        R_axis, Z_axis = get_magnetic_axis()
    finally:
        efit_to_boozer.deinit()
    assert np.isfinite(R_axis)
    assert np.isfinite(Z_axis)
    _debug("done test_get_magnetic_axis")

if __name__ == "__main__":
    pytest.main([__file__, "-s"])
