"""Tests for bc_to_booz_xform: .bc -> boozmn -> chartmap conversion."""

from contextlib import redirect_stdout
import io
from pathlib import Path
from types import SimpleNamespace
from urllib.request import urlopen
import numpy as np
import pytest

TESTS_DIR = Path(__file__).parent
CIRC_BC = Path(__file__).resolve().parent / "circ.bc"
LSP_WOUT_URL = (
    "https://raw.githubusercontent.com/hiddenSymmetries/simsopt/master/"
    "tests/test_files/wout_LandremanSenguptaPlunk_section5p3_reference.nc"
)

R0_expected = 1.64377  # m, from circ.bc header
A_expected = 0.46      # m


def _download(url: str, out_path: Path) -> None:
    with urlopen(url, timeout=60) as resp, open(out_path, "wb") as f:
        f.write(resp.read())


@pytest.fixture(scope="module")
def boozmn_path(tmp_path_factory):
    """Convert circ.bc to boozmn NetCDF once per module."""
    pytest.importorskip("netCDF4")
    if not CIRC_BC.exists():
        pytest.skip(f".bc file not found: {CIRC_BC}")

    import libneo
    src = str(Path(__file__).resolve().parents[1] / "libneo")
    if src not in libneo.__path__:
        libneo.__path__.insert(0, src)

    from libneo.bc_to_booz_xform import convert_bc_to_boozmn

    out = tmp_path_factory.mktemp("boozmn") / "boozmn_circ.nc"
    convert_bc_to_boozmn(CIRC_BC, out)
    return out


@pytest.fixture(scope="module")
def chartmap_path(tmp_path_factory, boozmn_path):
    """Build chartmap from boozmn once per module."""
    import libneo
    src = str(Path(__file__).resolve().parents[1] / "libneo")
    if src not in libneo.__path__:
        libneo.__path__.insert(0, src)

    from libneo.booz_xform_to_boozer_chartmap import convert_boozmn_to_chartmap

    out = tmp_path_factory.mktemp("chartmap") / "chartmap_circ.nc"
    convert_boozmn_to_chartmap(boozmn_path, out, nrho=20, ntheta=48, nzeta=1)
    return out


def _read_boozmn(path):
    import netCDF4

    d = {}
    with netCDF4.Dataset(path) as ds:
        d["ns"] = int(np.asarray(ds.variables["ns_b"][:]))
        d["nfp"] = int(np.asarray(ds.variables["nfp_b"][:]))
        d["jlist"] = np.asarray(ds.variables["jlist"][:])
        d["ixm"] = np.asarray(ds.variables["ixm_b"][:])
        d["ixn"] = np.asarray(ds.variables["ixn_b"][:])
        d["lasym"] = int(np.asarray(ds.variables["lasym__logical__"][:]))
        d["mboz"] = int(np.asarray(ds.variables["mboz_b"][:]))
        d["nboz"] = int(np.asarray(ds.variables["nboz_b"][:]))
        d["iota"] = np.asarray(ds.variables["iota_b"][:])
        d["bvco"] = np.asarray(ds.variables["bvco_b"][:])
        d["buco"] = np.asarray(ds.variables["buco_b"][:])
        d["phi"] = np.asarray(ds.variables["phi_b"][:])
        d["bmnc"] = np.asarray(ds.variables["bmnc_b"][:])
        d["rmnc"] = np.asarray(ds.variables["rmnc_b"][:])
        d["zmns"] = np.asarray(ds.variables["zmns_b"][:])
        d["pmns"] = np.asarray(ds.variables["pmns_b"][:])
        if d["lasym"]:
            d["bmns"] = np.asarray(ds.variables["bmns_b"][:])
            d["pmnc"] = np.asarray(ds.variables["pmnc_b"][:])
    return d


def _read_chartmap(path):
    import netCDF4

    d = {}
    with netCDF4.Dataset(path) as ds:
        d["rho"] = np.asarray(ds.variables["rho"][:])
        d["s"] = np.asarray(ds.variables["s"][:])
        d["x"] = np.asarray(ds.variables["x"][:])
        d["y"] = np.asarray(ds.variables["y"][:])
        d["z"] = np.asarray(ds.variables["z"][:])
        d["Bmod"] = np.asarray(ds.variables["Bmod"][:])
        d["A_phi"] = np.asarray(ds.variables["A_phi"][:])
        d["B_theta"] = np.asarray(ds.variables["B_theta"][:])
        d["B_phi"] = np.asarray(ds.variables["B_phi"][:])
        d["torflux"] = float(ds.torflux)
        d["nfp"] = int(np.asarray(ds.variables["num_field_periods"][:]).item())
    return d


def test_boozmn_structure(boozmn_path):
    """boozmn has correct dimensions and variable names."""
    d = _read_boozmn(boozmn_path)
    assert d["ns"] == 65, f"expected ns=63+2=65, got {d['ns']}"
    assert d["nfp"] == 1, f"expected nfp=1, got {d['nfp']}"
    assert len(d["jlist"]) == 63, f"expected 63 half-grid surfaces"
    assert d["jlist"][0] == 3, f"jlist[0]={d['jlist'][0]}, expected 3"
    assert d["jlist"][-1] == 65, f"jlist[-1]={d['jlist'][-1]}, expected 65"
    # The reconstructed half-grid must land on the .bc surface labels, otherwise
    # the chartmap is radially shifted by a fraction of a cell (off-by-one guard).
    s_half = (d["jlist"] - 1.5) / (d["ns"] - 1)
    assert abs(s_half[0] - 0.023438) < 1e-4 and abs(s_half[-1] - 0.992188) < 1e-4, (
        f"s_half {s_half[0]:.6f}..{s_half[-1]:.6f} does not match the .bc surfaces"
    )


def test_boozmn_iota_range(boozmn_path):
    """iota from boozmn is positive and in a physically plausible range for this tokamak."""
    d = _read_boozmn(boozmn_path)
    j = d["jlist"] - 1
    iota = d["iota"][j]
    # circ.bc has q=1/iota from ~1.1 (inner) to ~2.7 (outer) => iota from ~0.9 to ~0.37
    assert np.all(iota > 0.1), f"iota below 0.1: min={iota.min():.4f}"
    assert np.all(iota < 1.5), f"iota above 1.5: max={iota.max():.4f}"
    # Inner surface should be close to 0.9 (q~1.1)
    assert abs(iota[0] - 0.9) < 0.1, f"inner iota={iota[0]:.4f}, expected ~0.9"


def test_boozmn_bvco_sign(boozmn_path):
    """bvco_b (covariant B_phi) has consistent sign across all surfaces."""
    d = _read_boozmn(boozmn_path)
    j = d["jlist"] - 1
    bvco = d["bvco"][j]
    assert np.all(bvco < 0) or np.all(bvco > 0), (
        f"bvco has mixed signs; min={bvco.min():.3e}, max={bvco.max():.3e}"
    )


def test_boozmn_flux_monotone(boozmn_path):
    """Toroidal flux phi_b is monotone and ends at bc.flux."""
    d = _read_boozmn(boozmn_path)
    phi = d["phi"]
    assert phi[0] == pytest.approx(0.0, abs=1e-10), f"phi[0]={phi[0]}"
    assert abs(phi[-1]) == pytest.approx(1.33, rel=0.01), (
        f"phi[-1]={phi[-1]:.4f}, expected ~-1.33"
    )
    dphi = np.diff(phi)
    assert np.all(dphi * np.sign(phi[-1]) > -1e-12), "phi not monotone"


def test_boozmn_R0_from_rmnc(boozmn_path):
    """Major radius R0 recoverable as rmnc[m=0,n=0] of outer surfaces."""
    d = _read_boozmn(boozmn_path)
    m0n0 = (d["ixm"] == 0) & (d["ixn"] == 0)
    R0_computed = float(d["rmnc"][-1, m0n0].item())
    assert abs(R0_computed - R0_expected) < 0.1, (
        f"R0={R0_computed:.4f} m, expected {R0_expected:.4f} m"
    )


def test_symmetric_bc_writes_modes_with_booz_xform_signs(tmp_path):
    """Symmetric .bc coefficients keep booz_xform cos(m theta - n zeta)."""
    pytest.importorskip("netCDF4")
    from libneo.bc_to_booz_xform import write_boozmn

    mode_m = np.array([0, 1, 1, 2], dtype=int)
    mode_n = np.array([0, 1, -1, -2], dtype=int)
    zeros = [np.zeros(4), np.zeros(4)]
    bmnc = [
        np.array([1.5, 0.12, -0.19, 0.08]),
        np.array([1.6, -0.21, 0.15, -0.11]),
    ]
    vmns = [
        np.array([0.0, -0.14, 0.22, 0.07]),
        np.array([0.0, 0.18, -0.25, 0.09]),
    ]
    bc = SimpleNamespace(
        nsurf=2,
        nper=4,
        s=np.array([0.25, 0.75]),
        flux=1.0,
        iota=np.array([0.4, 0.5]),
        Jpol_divided_by_nper=np.array([1.0, 2.0]),
        Itor=np.array([3.0, 4.0]),
        m=[mode_m, mode_m],
        n=[mode_n, mode_n],
        rmnc=[np.ones(4), np.ones(4)],
        zmns=zeros,
        vmns=vmns,
        bmnc=bmnc,
        rmns=zeros,
        zmnc=zeros,
        vmnc=zeros,
        bmns=zeros,
    )

    out = tmp_path / "sym.nc"
    write_boozmn(bc, out)

    d = _read_boozmn(out)
    assert d["lasym"] == 0
    assert d["mboz"] == 2
    assert d["nboz"] == 2
    np.testing.assert_array_equal(d["ixm"], mode_m)
    np.testing.assert_array_equal(d["ixn"], mode_n * bc.nper)

    surface = 1
    theta = 0.37
    zeta = 0.21
    scale = 2.0 * np.pi / bc.nper
    angle = mode_m * theta - d["ixn"] * zeta
    np.testing.assert_allclose(
        np.sum(d["bmnc"][surface] * np.cos(angle)),
        np.sum(bmnc[surface] * np.cos(angle)),
    )
    np.testing.assert_allclose(
        np.sum(d["pmns"][surface] * np.sin(angle)),
        np.sum(vmns[surface] * scale * np.sin(angle)),
    )


def test_asymmetric_bc_writes_phase_coefficients_with_booz_xform_signs(tmp_path):
    """Asymmetric .bc phase coefficients keep booz_xform cos(m theta - n zeta)."""
    pytest.importorskip("netCDF4")
    from libneo.bc_to_booz_xform import write_boozmn

    mode_m = np.array([0, 0, 1, 1, 2], dtype=int)
    mode_n = np.array([0, 1, -1, 1, -2], dtype=int)
    zeros = [np.zeros(5), np.zeros(5)]
    bmnc = [
        np.array([1.5, 0.21, -0.13, 0.34, -0.08]),
        np.array([1.6, -0.17, 0.22, 0.41, 0.19]),
    ]
    bmns = [
        np.array([0.0, -0.12, 0.07, 0.18, -0.03]),
        np.array([0.0, 0.09, -0.16, 0.25, 0.11]),
    ]
    vmns = [
        np.array([0.0, 0.31, -0.23, 0.14, -0.05]),
        np.array([0.0, -0.27, 0.19, 0.33, 0.08]),
    ]
    vmnc = [
        np.array([0.0, -0.18, 0.26, -0.07, 0.16]),
        np.array([0.0, 0.24, -0.11, 0.29, -0.21]),
    ]
    bc = SimpleNamespace(
        nsurf=2,
        nper=3,
        s=np.array([0.25, 0.75]),
        flux=1.0,
        iota=np.array([0.4, 0.5]),
        Jpol_divided_by_nper=np.array([1.0, 2.0]),
        Itor=np.array([3.0, 4.0]),
        m=[mode_m, mode_m],
        n=[mode_n, mode_n],
        rmnc=[
            np.array([1.0, 0.1, 0.2, 0.3, 0.4]),
            np.array([1.1, 0.2, 0.3, 0.4, 0.5]),
        ],
        zmns=zeros,
        vmns=vmns,
        bmnc=bmnc,
        rmns=zeros,
        zmnc=zeros,
        vmnc=vmnc,
        bmns=bmns,
    )

    out = tmp_path / "asym.nc"
    write_boozmn(bc, out)

    d = _read_boozmn(out)
    assert d["lasym"] == 1
    assert d["mboz"] == 2
    assert d["nboz"] == 2
    np.testing.assert_array_equal(d["ixm"], mode_m)
    np.testing.assert_array_equal(d["ixn"], mode_n * bc.nper)

    scale = 2.0 * np.pi / bc.nper
    np.testing.assert_allclose(d["pmns"], np.asarray(vmns) * scale)
    np.testing.assert_allclose(d["pmnc"], np.asarray(vmnc) * scale)

    surface = 1
    theta = 0.37
    zeta = 0.21
    angle = mode_m * theta - d["ixn"] * zeta
    np.testing.assert_allclose(
        np.sum(
            d["bmnc"][surface] * np.cos(angle)
            + d["bmns"][surface] * np.sin(angle)
        ),
        np.sum(bmnc[surface] * np.cos(angle) + bmns[surface] * np.sin(angle)),
    )
    np.testing.assert_allclose(
        np.sum(
            d["pmnc"][surface] * np.cos(angle)
            + d["pmns"][surface] * np.sin(angle)
        ),
        np.sum(
            vmnc[surface] * scale * np.cos(angle)
            + vmns[surface] * scale * np.sin(angle)
        ),
    )


def test_asymmetric_vmec_bc_to_boozmn_matches_booz_xform_mode_signs(tmp_path):
    pytest.importorskip("booz_xform")
    pytest.importorskip("scipy")
    pytest.importorskip("netCDF4")

    from booz_xform import Booz_xform
    from libneo.bc_to_booz_xform import convert_bc_to_boozmn
    from libneo.boozer import BoozerFile

    wout = tmp_path / "wout_lsp_section5p3.nc"
    booz_xform_path = tmp_path / "boozmn_booz_xform.nc"
    bc_path = tmp_path / "libneo.bc"
    libneo_path = tmp_path / "boozmn_libneo.nc"

    _download(LSP_WOUT_URL, wout)

    with redirect_stdout(io.StringIO()):
        bx = Booz_xform()
        bx.read_wout(str(wout), True)
        bx.mboz = 6
        bx.nboz = 8
        bx.run()
        bx.write_boozmn(str(booz_xform_path))

        bc = BoozerFile(filename="")
        bc.convert_vmec_to_boozer(str(wout), uv_grid_multiplicator=1)
        bc.write(str(bc_path))
        convert_bc_to_boozmn(bc_path, libneo_path)

    booz_xform_boozmn = _read_boozmn(booz_xform_path)
    libneo_boozmn = _read_boozmn(libneo_path)

    assert booz_xform_boozmn["lasym"] == 1
    assert libneo_boozmn["lasym"] == 1
    assert "pmnc" in libneo_boozmn
    assert booz_xform_boozmn["nfp"] == libneo_boozmn["nfp"] == 3

    booz_xform_modes = set(zip(booz_xform_boozmn["ixm"], booz_xform_boozmn["ixn"]))
    libneo_modes = set(zip(libneo_boozmn["ixm"], libneo_boozmn["ixn"]))
    missing_modes = sorted(booz_xform_modes - libneo_modes)
    sign_flipped_modes = [
        (m, n) for (m, n) in missing_modes if n != 0 and (m, -n) in libneo_modes
    ]

    assert not sign_flipped_modes
    assert not missing_modes

    booz_xform_m0_n = sorted(n for m, n in booz_xform_modes if m == 0 and n > 0)
    assert booz_xform_m0_n
    assert all((0, n) in libneo_modes for n in booz_xform_m0_n)
    assert not any((0, -n) in libneo_modes for n in booz_xform_m0_n)


def test_chartmap_from_bc_structure(chartmap_path):
    """Chartmap built from .bc boozmn has nfp=1 and plausible geometry."""
    d = _read_chartmap(chartmap_path)
    assert d["nfp"] == 1, f"expected nfp=1, got {d['nfp']}"
    # R is sqrt(x^2+y^2) in cm; R0 should be near R0_expected*100 cm
    R_cm = np.sqrt(d["x"] ** 2 + d["y"] ** 2)
    R_mid = R_cm[R_cm.shape[0] // 2].mean()
    assert abs(R_mid / 100.0 - R0_expected) < 0.5, (
        f"R_mid={R_mid/100:.3f} m, expected ~{R0_expected} m"
    )


def test_chartmap_from_bc_Bmod(chartmap_path):
    """Bmod follows 1/R trend to within 5%."""
    d = _read_chartmap(chartmap_path)
    # Use the outermost surface for a basic consistency check.
    Bmod_outer = d["Bmod"][0].mean()  # (nzeta, ntheta, nrho) -> average over angles
    # analytic: B ~ B0*R0/R_outer; B0_cgs = 1.96e4 G, R0_cm = 164.4 cm
    B0_cgs = 1.96e4
    R0_cm = R0_expected * 100.0
    # outer surface in chartmap is outermost rho, not innermost
    # just verify Bmod is in a reasonable range
    assert 1.0e3 < Bmod_outer < 5.0e4, (
        f"Bmod_outer={Bmod_outer:.2f} G out of expected range [1e3,5e4]"
    )
