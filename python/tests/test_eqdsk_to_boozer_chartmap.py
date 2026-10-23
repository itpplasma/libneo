"""EFIT chartmaps from a consistent analytic Solov'ev equilibrium.

The historical EFIT31114 artifact has opposite signs in PSIRZ and the flux
header and is checked separately for diagnostic rejection.
"""

from pathlib import Path

import numpy as np
import pytest

from python.tests.efit_solovev import _write_solovev

TESTS_DIR = Path(__file__).resolve().parent
RESOURCES = TESTS_DIR.parents[1] / "test" / "resources"
EFIT = RESOURCES / "input_efit_file.dat"

METER_TO_CM = 1.0e2


def test_historical_flux_header_mismatch_is_rejected(tmp_path):
    from libneo.eqdsk_base import read_eqdsk
    from libneo.eqdsk_to_boozer_chartmap import convert_eqdsk_to_chartmap
    from scipy.interpolate import RectBivariateSpline

    eq = read_eqdsk(EFIT)
    psi = RectBivariateSpline(eq["R"], eq["Z"], eq["PsiVs"].T)
    # This public legacy artifact has grid-axis psi=+.108122192 while
    # SIMAG=-.108122148. Its LCFS similarly has the opposite sign to SIBRY.
    assert psi.ev(eq["Rpsi0"], eq["Zpsi0"]).item() > 0 > eq["PsiaxisVs"]
    assert np.median(psi.ev(eq["Lcfs"][:, 0], eq["Lcfs"][:, 1])) < 0 < eq["PsiedgeVs"]
    with pytest.raises(ValueError, match="not star-shaped|leaves the grid"):
        convert_eqdsk_to_chartmap(EFIT, tmp_path/"invalid.nc", nrho=8,
                                  ntheta=16, nzeta=8, ntheta_int=128, n_q=32)


@pytest.fixture(scope="module")
def chartmap_path(tmp_path_factory):
    pytest.importorskip("netCDF4")
    try:
        import _efit_to_boozer  # noqa: F401
        from efit_to_boozer.boozer import get_magnetic_axis  # noqa: F401
    except ImportError:
        pytest.skip("_efit_to_boozer not available")

    from libneo.eqdsk_to_boozer_chartmap import convert_eqdsk_to_chartmap

    out = tmp_path_factory.mktemp("chartmap") / "efit_boozer.nc"
    source = out.with_suffix(".geqdsk")
    _write_solovev(source)
    convert_eqdsk_to_chartmap(source, out)
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
    assert d["Bmod"].shape[2] >= 8, (
        f"expected a full toroidal grid, got {d['Bmod'].shape[2]} planes"
    )


def test_full_torus_geometry(chartmap_path):
    import netCDF4

    with netCDF4.Dataset(chartmap_path) as ds:
        rho = np.asarray(ds.variables["rho"][:])
        zeta = np.asarray(ds.variables["zeta"][:])
        x, y, z, b = [np.asarray(ds.variables[n][:]) for n in ["x", "y", "z", "Bmod"]]
    assert len(zeta) >= 8
    assert 0 <= rho[0] <= .01
    assert rho[-1] == 1
    xyz0 = (x[0]+1j*y[0])[None, ...]
    np.testing.assert_allclose(x+1j*y, xyz0*np.exp(1j*zeta[:, None, None]), atol=1e-9)
    np.testing.assert_allclose(z, np.broadcast_to(z[0], z.shape), atol=1e-9)
    np.testing.assert_allclose(b, np.broadcast_to(b[0], b.shape), atol=1e-9)


@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("toroidal_sign", [-1, 1])
def test_limiter_boundary_both_flux_directions(tmp_path, sign, toroidal_sign):
    from libneo.eqdsk_base import read_eqdsk, write_eqdsk
    from libneo.eqdsk_to_boozer_chartmap import convert_eqdsk_to_chartmap
    import netCDF4

    source = tmp_path/"solovev.geqdsk"
    _write_solovev(source)
    d = read_eqdsk(source)
    r = np.linspace(8.7, 11.3, 65)
    z = np.linspace(-1.3, 1.3, 65)
    theta = np.linspace(0, 2*np.pi, 129)
    d.update(nrgr=65, nzgr=65, rboxlength=2.6, zboxlength=2.6, R0=10.,
             rboxleft=8.7, zboxmid=0., Rpsi0=10., Zpsi0=0.,
             PsiaxisVs=0., PsiedgeVs=sign*.5, Btor_at_R0=2.*toroidal_sign, Ip=-sign*5e5,
             fprof=np.full(65, 20.*toroidal_sign), ptotprof=np.zeros(65),
             fdfdpsiprof=np.zeros(65), dpressdpsiprof=np.zeros(65),
             PsiVs=sign*.5*((r[None, :]-10)**2+z[:, None]**2),
             qprof=np.full(65, sign*toroidal_sign*2.), npbound=129, nplimiter=0,
             Lcfs=np.column_stack((10+np.cos(theta), np.sin(theta))),
             Limiter=np.empty((0, 2)))
    path = tmp_path/"circular.geqdsk"
    write_eqdsk(path, d)
    output = tmp_path/"chartmap.nc"
    convert_eqdsk_to_chartmap(path, output, nrho=24, ntheta=48, nstep=360,
        nlabel=32, ntheta_int=64, nsurfmax=400, nsurf=64, psimax=sign*.5e8,
        allow_poor_equilibrium=True, exp_bpol=0,
        override_reason="analytic circular geometry with intentionally zero GS profiles")
    with netCDF4.Dataset(output) as ds:
        x, y, z = [np.asarray(ds.variables[n][:]) for n in ["x", "y", "z"]]
        minor = np.sqrt((np.hypot(x, y)*.01-10)**2+(z*.01)**2)
        assert .98 < minor.max() < 1.01
        assert .98e8 < abs(ds.torflux) < 1.02e8
        # The exact Killing-coordinate covariant field is R*B_phi=20 T m.
        # magdata's normalization C also contains iota*I and is not G alone.
        np.testing.assert_allclose(ds.variables["B_phi"][:], 20e6*toroidal_sign, rtol=2e-5)
        assert np.sign(ds.torflux) == -toroidal_sign
        radius = np.hypot(x, y)*.01
        minor2 = (radius-10)**2+(z*.01)**2
        expected_b = np.sqrt((20/radius)**2+minor2/radius**2)
        np.testing.assert_allclose(ds.variables["Bmod"][:]*1e-4, expected_b, rtol=2e-5)


def test_bmod_positive(chartmap_path):
    d = _read_chartmap(chartmap_path)
    assert np.all(np.isfinite(d["Bmod"]))
    assert np.all(d["Bmod"] > 0.0), f"Bmod has non-positive values: min={d['Bmod'].min()}"


def test_geometry_follows_input_flux_surfaces(chartmap_path):
    from libneo.eqdsk_base import read_eqdsk
    from scipy.interpolate import RectBivariateSpline

    eq = read_eqdsk(chartmap_path.with_suffix(".geqdsk"))
    d = _read_chartmap(chartmap_path)
    r = np.hypot(d["x"][..., 0], d["y"][..., 0])*.01
    z = d["z"][..., 0]*.01
    flux = RectBivariateSpline(eq["R"], eq["Z"], eq["PsiVs"].T).ev(r, z)
    interior = (d["rho"] >= .1) & (d["rho"] <= .8)
    relative_spread = np.ptp(flux[interior], axis=1)/abs(eq["PsiedgeVs"]-eq["PsiaxisVs"])
    assert np.max(relative_spread) < 1e-3


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
    Ip = abs(float(read_eqdsk(chartmap_path.with_suffix(".geqdsk"))["Ip"]))
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
    # In geometric counter-clockwise theta, A_phi is the input COCOS-3
    # poloidal flux relative to the axis. Its sign is independent of B_tor.
    from libneo.eqdsk_base import read_eqdsk
    eq = read_eqdsk(chartmap_path.with_suffix(".geqdsk"))
    d = _read_chartmap(chartmap_path)
    A = d["A_phi"]
    assert np.sign(A[-1]) == np.sign(eq["PsiedgeVs"]-eq["PsiaxisVs"])
    assert abs(A[-1]) > abs(A[0]), (
        f"|A_phi| not increasing: |A_phi[0]|={abs(A[0]):.3e}, |A_phi[-1]|={abs(A[-1]):.3e}"
    )
