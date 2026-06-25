"""Cross-path field equality test: .bc direct Fourier evaluation vs chartmap.

Converts circ.bc -> boozmn -> chartmap, then compares |B| at a set of
(s, theta_B) points between:
  - direct Fourier summation from the .bc harmonics
  - the booz_xform_to_boozer_chartmap interpolation

This verifies that the .bc -> boozmn -> chartmap pipeline preserves field
values to numerical precision of the spline interpolation.
"""

from pathlib import Path
import numpy as np
import pytest
from scipy.interpolate import CubicSpline

CIRC_BC = Path(__file__).resolve().parent / "circ.bc"

TWOPI = 2.0 * np.pi


def _direct_bmod_from_bc(bc, s_test, theta_b_test):
    """Evaluate |B| directly from .bc Fourier harmonics at (s, theta_B) points.

    For n0b=0 (axisymmetric): bmnc[m,0] only; |B| = sum_m bmnc[m]*cos(m*theta_B).
    Spline in s between the .bc surfaces.
    """
    nsurf = bc.nsurf
    s_half = np.array(bc.s)
    m0 = np.array(bc.m[0], dtype=int)
    nmn = len(m0)

    # Build bmnc array (nsurf, nmn).
    bmnc = np.array([bc.bmnc[k] for k in range(nsurf)])

    n_test = len(s_test)
    bmod = np.zeros(n_test)
    for i in range(n_test):
        # Radial spline at s_test[i] for each mode.
        bmnc_at_s = np.array([
            float(CubicSpline(s_half, bmnc[:, j], extrapolate=True)(s_test[i]))
            for j in range(nmn)
        ])
        bmod[i] = float(np.dot(bmnc_at_s, np.cos(m0 * theta_b_test[i])))

    return bmod


@pytest.fixture(scope="module")
def chartmap_path(tmp_path_factory):
    """Build the .bc -> boozmn -> chartmap pipeline."""
    pytest.importorskip("netCDF4")
    if not CIRC_BC.exists():
        pytest.skip(f".bc file not found: {CIRC_BC}")

    import libneo
    src = str(Path(__file__).resolve().parents[1] / "libneo")
    if src not in libneo.__path__:
        libneo.__path__.insert(0, src)

    from libneo.bc_to_booz_xform import convert_bc_to_boozmn
    from libneo.booz_xform_to_boozer_chartmap import convert_boozmn_to_chartmap

    tmp = tmp_path_factory.mktemp("cross")
    boozmn = tmp / "boozmn_circ.nc"
    chartmap = tmp / "chartmap_circ.nc"
    convert_bc_to_boozmn(CIRC_BC, boozmn)
    convert_boozmn_to_chartmap(boozmn, chartmap, nrho=40, ntheta=96, nzeta=1)
    return chartmap


def test_bmod_cross_path(chartmap_path):
    """Bmod from .bc Fourier sum matches chartmap interpolation to 5e-3 relative.

    The error budget: the booz_xform_to_boozer_chartmap spline in rho introduces
    ~1-3e-3 error at mid-radius and up to 1% near the axis (s < 0.1) due to the
    near-axis power-law extrapolation for non-zero poloidal modes.  Test points
    are restricted to 0.2 < s < 0.9 where the error is bounded by 5e-3.
    """
    import netCDF4
    from scipy.interpolate import RegularGridInterpolator
    from libneo.boozer import BoozerFile

    bc = BoozerFile(str(CIRC_BC))

    with netCDF4.Dataset(chartmap_path) as ds:
        rho = np.asarray(ds.variables["rho"][:])
        theta = np.asarray(ds.variables["theta"][:])
        # Bmod has shape (nzeta, ntheta, nrho) in the NetCDF
        Bmod_nc = np.asarray(ds.variables["Bmod"][:])  # shape (nzeta, ntheta, nrho)
        # For nzeta=1: pick zeta=0 slice
        Bmod_2d = Bmod_nc[0, :, :].T  # shape (nrho, ntheta)

    # Build a 2D interpolator Bmod(rho, theta).
    theta_wrapped = np.append(theta, TWOPI)
    Bmod_wrapped = np.concatenate([Bmod_2d, Bmod_2d[:, :1]], axis=1)
    interp = RegularGridInterpolator(
        (rho, theta_wrapped),
        Bmod_wrapped,
        method="cubic",
        bounds_error=False,
        fill_value=None,
    )

    # Test points: 20 random (s, theta_B) in the mid-radius band.
    rng = np.random.default_rng(42)
    n_test = 20
    s_test = rng.uniform(0.20, 0.90, n_test)
    rho_test = np.sqrt(s_test)
    theta_b_test = rng.uniform(0.0, TWOPI, n_test)

    # Direct Fourier evaluation from .bc (in T).
    bmod_direct = _direct_bmod_from_bc(bc, s_test, theta_b_test)

    # Chartmap evaluation (in G): convert to T.
    bmod_chartmap_G = interp(np.column_stack([rho_test, theta_b_test]))
    bmod_chartmap = bmod_chartmap_G / 1.0e4

    for i in range(n_test):
        rel_err = abs(bmod_direct[i] - bmod_chartmap[i]) / abs(bmod_direct[i])
        assert rel_err < 5e-3, (
            f"point {i}: s={s_test[i]:.3f}, theta={theta_b_test[i]:.3f} rad; "
            f"|B|_direct={bmod_direct[i]:.5f} T, |B|_chartmap={bmod_chartmap[i]:.5f} T, "
            f"rel_err={rel_err:.2e}"
        )
