"""Tests that VMEC and chartmap R,Z evaluations are consistent.

When converting coordinates between VMEC and chartmap, the physical
position (R,Z) must be preserved exactly. This tests the coordinate
conversion functions used in the SQUID workflow.
"""
from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pytest


def _download(url: str, out_path: Path) -> None:
    with urlopen(url, timeout=60) as resp, open(out_path, "wb") as f:
        f.write(resp.read())


def _vmec_coords_s(
    wout_path: Path, s: np.ndarray, theta: np.ndarray, zeta: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Compute VMEC (R,Z) from (s, theta, zeta)."""
    from netCDF4 import Dataset

    with Dataset(wout_path, "r") as ds:
        rmnc = np.array(ds.variables["rmnc"][:], dtype=float)
        zmns = np.array(ds.variables["zmns"][:], dtype=float)
        xm = np.array(ds.variables["xm"][:], dtype=float)
        xn = np.array(ds.variables["xn"][:], dtype=float)
        ns = int(rmnc.shape[0])

    s_clamped = np.clip(s.astype(float), 0.0, 1.0)
    u = s_clamped * float(ns - 1)
    js0 = np.floor(u).astype(int)
    js0 = np.clip(js0, 0, ns - 1)
    js1 = np.clip(js0 + 1, 0, ns - 1)
    t = u - js0.astype(float)

    angle = theta[:, None] * xm[None, :] - zeta[:, None] * xn[None, :]
    c = np.cos(angle)
    sng = np.sin(angle)

    rm = (1.0 - t)[:, None] * rmnc[js0, :] + t[:, None] * rmnc[js1, :]
    zm = (1.0 - t)[:, None] * zmns[js0, :] + t[:, None] * zmns[js1, :]
    R = np.sum(rm * c, axis=1)
    Z = np.sum(zm * sng, axis=1)
    return R, Z


def _chartmap_rz(
    chartmap_path: Path,
    rho: np.ndarray,
    theta: np.ndarray,
    zeta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute chartmap (R,Z) from (rho, theta, zeta) via trilinear interpolation."""
    from netCDF4 import Dataset

    with Dataset(chartmap_path, "r") as ds:
        rho_grid = np.array(ds.variables["rho"][:], dtype=float)
        theta_grid = np.array(ds.variables["theta"][:], dtype=float)
        zeta_grid = np.array(ds.variables["zeta"][:], dtype=float)
        x = np.array(ds.variables["x"][:, :, :], dtype=float)
        y = np.array(ds.variables["y"][:, :, :], dtype=float)
        z = np.array(ds.variables["z"][:, :, :], dtype=float)
        nfp = int(np.array(ds.variables["num_field_periods"][...]))
        units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
        scale = 0.01 if units == "cm" else 1.0

    r_grid = np.sqrt((x * scale) ** 2 + (y * scale) ** 2)
    z_grid = z * scale

    rho = np.clip(rho.astype(float), float(rho_grid[0]), float(rho_grid[-1]))
    theta = np.mod(theta.astype(float), 2.0 * np.pi)
    zeta_period = 2.0 * np.pi / float(nfp)
    zeta = np.mod(zeta.astype(float), zeta_period)

    dr = float(rho_grid[1] - rho_grid[0])
    dth = float(2.0 * np.pi) / float(theta_grid.size)
    dz = float(zeta_period) / float(zeta_grid.size)

    ur = (rho - float(rho_grid[0])) / dr
    ut = theta / dth
    uz = zeta / dz

    ir0 = np.floor(ur).astype(int)
    it0 = np.floor(ut).astype(int) % theta_grid.size
    iz0 = np.floor(uz).astype(int) % zeta_grid.size

    ir0 = np.clip(ir0, 0, rho_grid.size - 2)
    ir1 = ir0 + 1
    it1 = (it0 + 1) % theta_grid.size
    iz1 = (iz0 + 1) % zeta_grid.size

    tr = ur - ir0.astype(float)
    tt = ut - np.floor(ut)
    tz = uz - np.floor(uz)

    def tri(arr):
        c000 = arr[iz0, it0, ir0]
        c001 = arr[iz0, it0, ir1]
        c010 = arr[iz0, it1, ir0]
        c011 = arr[iz0, it1, ir1]
        c100 = arr[iz1, it0, ir0]
        c101 = arr[iz1, it0, ir1]
        c110 = arr[iz1, it1, ir0]
        c111 = arr[iz1, it1, ir1]
        v00 = c000 * (1.0 - tr) + c001 * tr
        v01 = c010 * (1.0 - tr) + c011 * tr
        v10 = c100 * (1.0 - tr) + c101 * tr
        v11 = c110 * (1.0 - tr) + c111 * tr
        v0 = v00 * (1.0 - tt) + v01 * tt
        v1 = v10 * (1.0 - tt) + v11 * tt
        return v0 * (1.0 - tz) + v1 * tz

    R = tri(r_grid)
    Z = tri(z_grid)
    return R, Z


def _invert_chartmap(
    chartmap_path: Path,
    R_target: np.ndarray,
    Z_target: np.ndarray,
    zeta: np.ndarray,
    rho_init: np.ndarray,
    theta_init: np.ndarray,
    max_iter: int = 30,
    tol_m: float = 1.0e-8,
) -> tuple[np.ndarray, np.ndarray]:
    """Invert chartmap to find (rho, theta) for given (R, Z, zeta)."""
    n = len(R_target)
    rho = np.clip(rho_init.copy(), 0.0, 0.999)
    theta = theta_init.copy() % (2.0 * np.pi)

    for _ in range(max_iter):
        R, Z = _chartmap_rz(chartmap_path, rho, theta, zeta)
        f1 = R - R_target
        f2 = Z - Z_target
        err = np.maximum(np.abs(f1), np.abs(f2))
        if np.all(err < tol_m):
            break

        dr = 1.0e-4
        dth = 1.0e-4
        Rr, Zr = _chartmap_rz(chartmap_path, np.clip(rho + dr, 0, 0.999), theta, zeta)
        Rt, Zt = _chartmap_rz(chartmap_path, rho, theta + dth, zeta)

        dR_dr = (Rr - R) / dr
        dZ_dr = (Zr - Z) / dr
        dR_dth = (Rt - R) / dth
        dZ_dth = (Zt - Z) / dth

        det = dR_dr * dZ_dth - dR_dth * dZ_dr
        det = np.where(np.abs(det) < 1e-16, 1e-16, det)

        drho = (dZ_dth * f1 - dR_dth * f2) / det
        dtheta = (-dZ_dr * f1 + dR_dr * f2) / det

        rho = np.clip(rho - 0.5 * drho, 0.0, 0.999)
        theta = (theta - 0.5 * dtheta) % (2.0 * np.pi)

    R_final, Z_final = _chartmap_rz(chartmap_path, rho, theta, zeta)
    err_R = np.abs(R_final - R_target)
    err_Z = np.abs(Z_final - Z_target)
    max_err = float(np.max(np.maximum(err_R, err_Z)))
    if max_err >= tol_m:
        raise RuntimeError(
            f"chartmap inversion did not converge within {max_iter} iterations. "
            f"Max error: {max_err:.2e} m"
        )

    return rho, theta


@pytest.mark.network
def test_vmec_to_chartmap_rz_roundtrip(tmp_path: Path) -> None:
    """Test that VMEC -> R,Z -> chartmap inversion preserves R,Z."""
    pytest.importorskip("map2disc")
    pytest.importorskip("shapely")

    from libneo.chartmap import write_chartmap_from_vmec_boundary

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    chartmap = tmp_path / "chartmap.nc"
    write_chartmap_from_vmec_boundary(
        wout,
        chartmap,
        nrho=33,
        ntheta=65,
        nzeta=33,
        boundary_offset=0.0,
        boundary_param="arc",
    )

    # Generate test points in VMEC coordinates
    n_points = 20
    np.random.seed(42)
    s_vmec = 0.2 + 0.6 * np.random.rand(n_points)  # s in [0.2, 0.8]
    theta_vmec = 2 * np.pi * np.random.rand(n_points)
    zeta_vmec = np.pi / 6 * np.random.rand(n_points)  # within one field period

    # Compute R,Z from VMEC
    R_vmec, Z_vmec = _vmec_coords_s(wout, s_vmec, theta_vmec, zeta_vmec)

    # Invert to find chartmap coordinates
    rho_init = np.sqrt(s_vmec)
    theta_init = theta_vmec
    rho_cm, theta_cm = _invert_chartmap(
        chartmap, R_vmec, Z_vmec, zeta_vmec, rho_init, theta_init, max_iter=50
    )

    # Compute R,Z from chartmap
    R_cm, Z_cm = _chartmap_rz(chartmap, rho_cm, theta_cm, zeta_vmec)

    # R,Z should match
    err_R = np.abs(R_vmec - R_cm)
    err_Z = np.abs(Z_vmec - Z_cm)
    max_err = max(np.max(err_R), np.max(err_Z))

    assert max_err < 1e-4, f"R,Z roundtrip error {max_err:.2e} m exceeds tolerance"


@pytest.mark.network
def test_chartmap_theta_param_preserves_vmec_theta(tmp_path: Path) -> None:
    """With theta parameterization, chartmap theta should equal VMEC theta."""
    pytest.importorskip("map2disc")
    pytest.importorskip("shapely")

    from libneo.chartmap import write_chartmap_from_vmec_boundary

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    chartmap_theta = tmp_path / "chartmap_theta.nc"
    write_chartmap_from_vmec_boundary(
        wout,
        chartmap_theta,
        nrho=33,
        ntheta=65,
        nzeta=33,
        boundary_offset=0.0,
        boundary_param="theta",
    )

    # Test points
    n_points = 20
    np.random.seed(123)
    s_vmec = 0.3 + 0.5 * np.random.rand(n_points)
    theta_vmec = 2 * np.pi * np.random.rand(n_points)
    zeta_vmec = np.pi / 6 * np.random.rand(n_points)

    R_vmec, Z_vmec = _vmec_coords_s(wout, s_vmec, theta_vmec, zeta_vmec)

    rho_init = np.sqrt(s_vmec)
    rho_cm, theta_cm = _invert_chartmap(
        chartmap_theta, R_vmec, Z_vmec, zeta_vmec, rho_init, theta_vmec, max_iter=50
    )

    # With theta parameterization, the inverted theta should match VMEC theta
    # (modulo numerical tolerance from inversion)
    theta_diff = np.abs(theta_cm - (theta_vmec % (2 * np.pi)))
    theta_diff = np.minimum(theta_diff, 2 * np.pi - theta_diff)  # handle wrap

    # Allow some tolerance due to grid interpolation and Newton inversion
    # With 65 theta points, expect ~0.1 rad grid spacing, so 0.2 rad tolerance is reasonable
    assert np.max(theta_diff) < 0.2, (
        f"Theta mismatch {np.max(theta_diff):.3f} rad with theta parameterization"
    )


@pytest.mark.network
def test_chartmap_arc_param_changes_theta(tmp_path: Path) -> None:
    """With arc parameterization, chartmap theta differs from VMEC theta."""
    pytest.importorskip("map2disc")
    pytest.importorskip("shapely")

    from libneo.chartmap import write_chartmap_from_vmec_boundary

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    chartmap_arc = tmp_path / "chartmap_arc.nc"
    write_chartmap_from_vmec_boundary(
        wout,
        chartmap_arc,
        nrho=33,
        ntheta=65,
        nzeta=33,
        boundary_offset=0.0,
        boundary_param="arc",
    )

    # Test points
    n_points = 20
    np.random.seed(456)
    s_vmec = 0.3 + 0.5 * np.random.rand(n_points)
    theta_vmec = 2 * np.pi * np.random.rand(n_points)
    zeta_vmec = np.pi / 6 * np.random.rand(n_points)

    R_vmec, Z_vmec = _vmec_coords_s(wout, s_vmec, theta_vmec, zeta_vmec)

    rho_init = np.sqrt(s_vmec)
    rho_cm, theta_cm = _invert_chartmap(
        chartmap_arc, R_vmec, Z_vmec, zeta_vmec, rho_init, theta_vmec, max_iter=50
    )

    # With arc parameterization, R,Z must still match but theta will differ
    R_cm, Z_cm = _chartmap_rz(chartmap_arc, rho_cm, theta_cm, zeta_vmec)
    max_rz_err = max(np.max(np.abs(R_vmec - R_cm)), np.max(np.abs(Z_vmec - Z_cm)))
    assert max_rz_err < 1e-4, f"R,Z error {max_rz_err:.2e} m with arc parameterization"

    # Theta values should be different (arclength reparameterization)
    theta_diff = np.abs(theta_cm - (theta_vmec % (2 * np.pi)))
    theta_diff = np.minimum(theta_diff, 2 * np.pi - theta_diff)

    # Expect significant theta differences for non-circular cross-sections
    # NCSX has shaped cross-sections so theta should differ noticeably
    mean_diff = np.mean(theta_diff)
    assert mean_diff > 0.01, (
        f"Expected theta difference with arc param, got mean diff {mean_diff:.4f} rad"
    )
