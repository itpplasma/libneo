"""Strict EFIT-to-BC export checked against synthetic analytic field integrals."""
import importlib
import importlib.util
import os
from pathlib import Path
import sys

import numpy as np
import pytest

from python.tests.efit_solovev import (
    A, B, MU0, PSI_EDGE, _axis, _psi, _grad,
    _contour_oracle, _write_solovev,
)


NSURF, MPOL, NTHETA, NQ = 6, 20, 128, 32
S = (np.arange(NSURF)+.5)/NSURF
SPECTRA = ("rmnc", "rmns", "zmnc", "zmns", "vmnc", "vmns", "bmnc", "bmns")


@pytest.fixture
def export_api(monkeypatch):
    """Use the requested compiled revision and this checkout's Python front end."""
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    monkeypatch.setenv("OPENBLAS_NUM_THREADS", "1")
    monkeypatch.setenv("MPLBACKEND", "Agg")
    root = Path(__file__).resolve().parents[2]
    monkeypatch.syspath_prepend(str(root/"python"))
    extension_path = os.environ.get("LIBNEO_EFIT_EXTENSION")
    if extension_path:
        spec = importlib.util.spec_from_file_location("_efit_to_boozer", extension_path)
        extension = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(extension)
    else:
        extension = pytest.importorskip("_efit_to_boozer")
    monkeypatch.setitem(sys.modules, "_efit_to_boozer", extension)
    spec = importlib.util.spec_from_file_location(
        "efit_to_boozer.boozer", root/"python/efit_to_boozer/boozer.py")
    contour = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(contour)
    monkeypatch.setitem(sys.modules, "efit_to_boozer.boozer", contour)
    exporter = importlib.import_module("libneo.eqdsk_to_boozer")
    return exporter.convert_eqdsk_to_bc, contour, extension.efit_to_boozer


@pytest.fixture(scope="module")
def toroidal_flux():
    from scipy.integrate import quad

    # Positive flux per radian, independent of either input field polarity.
    return quad(lambda p: _contour_oracle(p)[0], _psi(*_axis()), PSI_EDGE)[0]


def _polar_volume_density(psi):
    """Return <R rho / psi_rho> and RMS rho using only the analytic polynomial."""
    angle = np.arange(1024)*2*np.pi/1024
    co, si = np.cos(angle), np.sin(angle)
    ra, za = _axis()
    lo, hi = np.zeros(angle.size), np.full(angle.size, 1.5)
    for _ in range(65):
        rho = (lo+hi)/2
        outside = _psi(ra+rho*co, za+rho*si) > psi
        hi, lo = np.where(outside, rho, hi), np.where(outside, lo, rho)
    rho = (lo+hi)/2
    r, z = ra+rho*co, za+rho*si
    pr, pz = _grad(r, z)
    return np.mean(r*rho/(pr*co+pz*si)), np.sqrt(np.mean(rho*rho))


def _fourier(bc, surface):
    theta = np.arange(128)*2*np.pi/128
    modes = np.asarray(bc.m[surface])[:, None]
    co, si = np.cos(modes*theta), np.sin(modes*theta)
    values, derivatives = [], []
    for cosine, sine in zip(SPECTRA[::2], SPECTRA[1::2]):
        c, s = np.asarray(getattr(bc, cosine)[surface]), np.asarray(getattr(bc, sine)[surface])
        values.append(c@co+s@si)
        derivatives.append(-(c*modes[:, 0])@si+(s*modes[:, 0])@co)
    return values, derivatives


@pytest.mark.parametrize("psi_sign,f_sign", [(1, 1), (-1, 1), (1, -1), (-1, -1)])
@pytest.mark.parametrize("route", [1, 2])
def test_bc_matches_analytic_equilibrium(tmp_path, export_api, toroidal_flux,
                                       psi_sign, f_sign, route):
    from libneo.boozer import BoozerFile

    convert, contour, _ = export_api
    gfile, output = tmp_path/"solovev.g", tmp_path/"equilibrium.bc"
    _write_solovev(gfile, psi_sign, f_sign)
    assert convert(gfile, output, nsurf=NSURF, mpol=MPOL, n_theta=NTHETA,
                   n_q=NQ, route=route) == output
    bc = BoozerFile(str(output))
    np.testing.assert_allclose(bc.s, S, rtol=0, atol=1e-16)
    assert bc.nper == 1 and bc.n0b == 0
    assert bc.flux == pytest.approx(f_sign*2*np.pi*toroidal_flux, rel=3e-6)
    assert bc.R == pytest.approx(_axis()[0], abs=2e-8)
    assert bc.a == pytest.approx(_polar_volume_density(PSI_EDGE)[1], rel=3e-7)
    metadata = contour.get_contour_transform_extended(
        gfile, S, MPOL, n_theta=NTHETA, n_q=NQ, route=route)
    # These comparisons expose precision loss in serialization independently
    # of the looser interpolation/truncation tolerances of the field oracle.
    scale = np.array([.01, .01, .01, .01, 1., 1., 1e-4, 1e-4])[:, None, None]
    saved = np.stack([np.asarray(getattr(bc, name)).T for name in SPECTRA])
    np.testing.assert_allclose(saved, metadata["spectra"]*scale, rtol=2e-15, atol=2e-16)
    np.testing.assert_allclose(bc.iota, metadata["iota"], rtol=2e-15)
    np.testing.assert_allclose(bc.Jpol_divided_by_nper, -5*metadata["curr_g"], rtol=2e-15)
    np.testing.assert_allclose(bc.Itor, -5*metadata["curr_i"], rtol=2e-15)
    np.testing.assert_allclose(bc.pprime, .1*metadata["pp_file"], rtol=2e-15)
    native_volume = -np.abs(metadata["sqrtg00"]*metadata["psitor_max"])*1e-6*(2*np.pi)**2
    np.testing.assert_allclose(bc.sqrt_g_00, native_volume, rtol=2e-15)
    for k in range(NSURF):
        (r, z, v, bmod), (rt, zt, vt, _) = _fourier(bc, k)
        psi = _psi(r, z)
        assert np.ptp(psi)/(PSI_EDGE-_psi(*_axis())) < 2e-7
        q, current_i, current_g = _contour_oracle(float(np.mean(psi)), psi_sign, f_sign)
        assert bc.iota[k] == pytest.approx(1/q, rel=3e-7)
        assert bc.Jpol_divided_by_nper[k] == pytest.approx(-5e6*current_g, rel=3e-7)
        assert bc.Itor[k] == pytest.approx(-5e6*current_i, rel=3e-7)
        pr, pz = _grad(r, z)
        np.testing.assert_allclose(bmod, np.sqrt(pr*pr+pz*pz+current_g**2)/r, rtol=3e-7)
        # Along a straight field line phi_B=theta_B/iota and
        # phi_cyl=phi_B+2*pi*v(theta_B), so the physical field is tangent.
        br, bz, bphi = -psi_sign*pz/r, psi_sign*pr/r, current_g/r
        tangent_scale = (br*rt+bz*zt)/(rt*rt+zt*zt)
        np.testing.assert_allclose(bphi, tangent_scale*r*(1/bc.iota[k]+2*np.pi*vt),
                                   rtol=3e-6, atol=1e-8)
        dpsi_ds = toroidal_flux/abs(q)
        pressure_derivative = -(8*B+2*A)/MU0*dpsi_ds
        assert bc.pprime[k] == pytest.approx(pressure_derivative, rel=3e-6)
        density, _ = _polar_volume_density(float(np.mean(psi)))
        independent_volume = -(2*np.pi)**2*density*dpsi_ds
        assert bc.sqrt_g_00[k] == pytest.approx(independent_volume, rel=3e-6)


def test_extended_contours_preserve_legacy_raw_outputs(tmp_path, export_api):
    _, contour, native = export_api
    gfile = tmp_path/"solovev.g"
    _write_solovev(gfile, psi_sign=-1)
    args = (str(gfile), S, MPOL, NTHETA, NQ, 1, False, 0., False, "")
    original = native.contour_transform(*args)
    extended = native.contour_transform_extended(*args)
    assert len(original) == 6 and len(extended) == 8
    for old, new in zip(original[:3], extended[:3]):
        np.testing.assert_array_equal(old, new)
    assert original[3:] == extended[5:]
    old = contour.get_contour_transform(gfile, S, MPOL, n_theta=NTHETA, n_q=NQ)
    new = contour.get_contour_transform_extended(gfile, S, MPOL, n_theta=NTHETA, n_q=NQ)
    assert set(new) == set(old) | {"pp_file", "sqrtg00", "rminor", "sig_theta"}
    for name, value in old.items():
        np.testing.assert_equal(value, new[name])
    np.testing.assert_array_equal(new["pp_file"], extended[3][0])
    np.testing.assert_array_equal(new["sqrtg00"], extended[3][1])
    assert new["rminor"] == extended[4]
    assert new["sig_theta"] == extended[2][3]


@pytest.mark.parametrize("failure", ["poor", "malformed", "nan_bound", "nan_grid", 3, 4, 5])
def test_rejected_export_preserves_existing_output(tmp_path, monkeypatch, export_api, failure):
    from libneo.eqdsk_base import read_eqdsk, write_eqdsk

    monkeypatch.chdir(tmp_path)
    convert, _, _ = export_api
    gfile, output = tmp_path/"input.g", tmp_path/"existing.bc"
    _write_solovev(gfile, pprime_scale=.7 if failure == "poor" else 1)
    options = dict(nsurf=NSURF, mpol=MPOL, n_theta=NTHETA, n_q=NQ)
    if failure == "malformed":
        gfile.write_text("truncated G-EQDSK\n")
    elif failure == "nan_bound":
        options["psimax"] = np.nan
    elif failure == "nan_grid":
        equilibrium = read_eqdsk(gfile)
        equilibrium["PsiVs"][32, 32] = np.nan
        write_eqdsk(gfile, equilibrium)
    elif isinstance(failure, int):
        options["route"] = failure
    original = b"previous output must survive a rejected conversion\n"
    output.write_bytes(original)
    before = set(tmp_path.rglob("*"))
    with pytest.raises(ValueError):
        convert(gfile, output, **options)
    assert output.read_bytes() == original
    assert set(tmp_path.rglob("*")) == before


def test_default_boozer_writers_preserve_historical_bytes(tmp_path, monkeypatch, export_api):
    from scipy.constants import mu_0

    writer = importlib.import_module("libneo.boozer")
    monkeypatch.setattr("getpass.getuser", lambda: "writer-test")
    output = tmp_path/"historical.bc"
    writer.write_boozer_head(output, "", 0, 1, 0, 1, 1, 2., 3., 4.)
    writer.append_boozer_block_head(output, .25, -.5, -3*mu_0/(2*np.pi),
                                   -4*mu_0/(2*np.pi), 5., -6/(4*np.pi**2), 1)
    writer.append_boozer_block(
        output, [0, 1], np.array([0, 0]),
        np.array([1-2j, -1+2j]), np.array([3-4j, -3+4j]),
        np.array([5-6j, -5+6j]), np.array([7-8j, -7+8j]), 1)
    expected = (
        b"CC Boozer-coordinate data file\n"
        b"CC Version: \n"
        b"CC Author: writer-test\n"
        b"CC shot:    0\n"
        b" m0b   n0b  nsurf  nper    flux [Tm^2]        a [m]          R [m]\n"
        b"   1     0    1    1  2.000000e+00   3.000000e+00   4.000000e+00\n"
        b"        s               iota           Jpol/nper          Itor            pprime         sqrt g(0,0)\n"
        b"                                          [A]           [A]             [Pa]         (dV/ds)/nper\n"
        b"   2.50000000e-01  -5.00000000e-01   3.00000000e+00   4.00000000e+00   5.00000000e+00   6.00000000e+00\n"
        b"    m    n      rmnc [m]         rmns [m]         zmnc [m]         zmns [m]         vmnc [ ]         vmns [ ]         bmnc [T]         bmns [T]\n"
        b"    0    0   1.00000000e+00   2.00000000e+00   3.00000000e+00   4.00000000e+00   5.00000000e+00   6.00000000e+00   7.00000000e+00   8.00000000e+00\n"
        b"    1    0  -1.00000000e+00  -2.00000000e+00  -3.00000000e+00  -4.00000000e+00  -5.00000000e+00  -6.00000000e+00  -7.00000000e+00  -8.00000000e+00\n"
    )
    assert output.read_bytes() == expected
    rewritten = tmp_path/"rewritten.bc"
    writer.BoozerFile(str(output)).write(rewritten)
    assert rewritten.read_bytes() == expected
