"""Tests for eqdsk_to_boozer_chartmap on the standard diverted equilibrium.

Drives efit_to_boozer with its standard settings (matching the committed
efit_to_boozer.inp) on test/resources/input_efit_file.dat, the same g-file the
efit_to_boozer test suite uses, and checks that the assembled Boozer chartmap is
structurally and physically consistent.
"""

from pathlib import Path

import numpy as np
import pytest

TESTS_DIR = Path(__file__).resolve().parent
RESOURCES = TESTS_DIR.parents[1] / "test" / "resources"
EFIT = RESOURCES / "input_efit_file.dat"
CONVEX = RESOURCES / "convexwall.dat"

METER_TO_CM = 1.0e2


@pytest.fixture(scope="module")
def chartmap_path(tmp_path_factory):
    pytest.importorskip("netCDF4")
    if not EFIT.exists():
        pytest.skip(f"equilibrium not found: {EFIT}")
    try:
        import _efit_to_boozer  # noqa: F401
        from efit_to_boozer.boozer import get_magnetic_axis  # noqa: F401
    except ImportError:
        pytest.skip("_efit_to_boozer not available")

    from libneo.eqdsk_to_boozer_chartmap import convert_eqdsk_to_chartmap

    out = tmp_path_factory.mktemp("chartmap") / "efit_boozer.nc"
    convert_eqdsk_to_chartmap(EFIT, out, convexwall=str(CONVEX))
    return out


def _read_chartmap(path):
    import netCDF4

    d = {}
    with netCDF4.Dataset(path) as ds:
        for k in ("rho", "s", "theta", "zeta", "A_phi", "B_theta", "B_phi"):
            d[k] = np.asarray(ds.variables[k][:])
        for k in ("Bmod", "x", "y", "z"):
            d[k] = np.transpose(np.asarray(ds.variables[k][:]), (2, 1, 0))
        d["torflux"] = float(ds.torflux)
        d["nfp"] = int(np.asarray(ds.variables["num_field_periods"][:]).item())
    return d


def test_axisymmetric_structure(chartmap_path):
    d = _read_chartmap(chartmap_path)
    assert d["nfp"] == 1, f"expected nfp=1, got {d['nfp']}"
    assert d["Bmod"].shape[2] == 1, (
        f"expected nzeta=1 for axisymmetric chartmap, got {d['Bmod'].shape[2]}"
    )


def test_bmod_positive(chartmap_path):
    d = _read_chartmap(chartmap_path)
    assert np.all(np.isfinite(d["Bmod"]))
    assert np.all(d["Bmod"] > 0.0), f"Bmod has non-positive values: min={d['Bmod'].min()}"


def test_geometry_plausible(chartmap_path):
    d = _read_chartmap(chartmap_path)
    R_m = np.sqrt(d["x"][..., 0] ** 2 + d["y"][..., 0] ** 2) / METER_TO_CM
    assert np.all(R_m > 0.0)
    R_axis = float(R_m[0].mean())
    assert 1.0 < R_axis < 2.5, f"axis R = {R_axis:.3f} m, expected a tokamak-scale value"


def test_covariant_components_tokamak_ordering(chartmap_path):
    # B_phi = G(s) = R*B_tor; B_theta = I(s) from the enclosed toroidal
    # current. For a tokamak I << G (NOT I = G/q: iota is the ratio of the
    # contravariant components, while I/G measures enclosed current).
    # Signs of both track the equilibrium's field and current direction.
    d = _read_chartmap(chartmap_path)
    assert np.all(np.isfinite(d["B_theta"])) and np.all(d["B_theta"] != 0.0), (
        f"B_theta (I) has zero/non-finite values: min|I|={np.abs(d['B_theta']).min():.4e}"
    )
    assert np.all(np.abs(d["B_phi"]) > 0.0), (
        f"B_phi (G) has zero values: min|G|={np.abs(d['B_phi']).min():.4e}"
    )
    ratio = np.abs(d["B_theta"]).max() / np.abs(d["B_phi"]).min()
    assert ratio < 0.5, (
        f"|I/G| = {ratio:.3f}: covariant B_theta is not small against B_phi, "
        "which suggests the fabricated I = G/q instead of the enclosed-current I"
    )


def test_B_theta_matches_plasma_current(chartmap_path):
    # Ampere: 2*pi*I(s) = mu0 * I_tor(s), so at the outermost stored surface
    # I approaches mu0*Ip/(2*pi). In CGS G*cm that is 0.2*Ip[A]. The stored
    # grid stops just inside the LCFS, so allow a one-sided tolerance.
    from libneo.eqdsk_base import read_eqdsk

    d = _read_chartmap(chartmap_path)
    Ip = abs(float(read_eqdsk(EFIT)["Ip"]))
    I_edge_expected = 0.2 * Ip  # G*cm
    I_edge = abs(float(d["B_theta"][-1]))
    assert 0.7 * I_edge_expected < I_edge <= 1.05 * I_edge_expected, (
        f"I at outermost surface = {I_edge:.4e} G*cm, expected close to "
        f"mu0*Ip/2pi = {I_edge_expected:.4e} G*cm"
    )


def test_torflux_nonzero(chartmap_path):
    d = _read_chartmap(chartmap_path)
    assert np.isfinite(d["torflux"]) and d["torflux"] != 0.0, (
        f"torflux={d['torflux']:.3e} must be finite and nonzero"
    )


def test_write_inp_psimax(tmp_path):
    # psimax bounds the flux-surface scan; required for equilibria without an
    # X-point, whose psi keeps rising to the box edge. Default stays 1e10.
    from libneo.eqdsk_to_boozer_chartmap import _write_inp

    inp = tmp_path / "efit_to_boozer.inp"
    _write_inp(inp, "dummy.geqdsk", psimax=4.25e8)
    assert "425000000.0" in inp.read_text()

    inp_default = tmp_path / "efit_to_boozer_default.inp"
    _write_inp(inp_default, "dummy.geqdsk")
    assert "10000000000.0" in inp_default.read_text()


def test_A_phi_consistent_with_torflux(chartmap_path):
    # A_phi(s) = -torflux * integral_0^s iota ds' with iota > 0, so sign(A_phi)
    # is opposite to sign(torflux) and |A_phi| increases with s. The absolute
    # signs track the equilibrium's field direction, not a fixed convention.
    d = _read_chartmap(chartmap_path)
    A = d["A_phi"]
    assert np.sign(A[-1]) == -np.sign(d["torflux"]), (
        f"sign(A_phi)={np.sign(A[-1])} inconsistent with -sign(torflux)"
    )
    assert abs(A[-1]) > abs(A[0]), (
        f"|A_phi| not increasing: |A_phi[0]|={abs(A[0]):.3e}, |A_phi[-1]|={abs(A[-1]):.3e}"
    )
