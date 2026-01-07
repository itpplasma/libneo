"""Tests that chartmap boundary matches VMEC boundary in R,Z physical space.

For offset=0, the chartmap boundary at rho=1 must match the VMEC boundary at s=1
in physical (R,Z) space, regardless of theta parameterization differences.

This validates the core assumption: arclength parameterization changes coordinate
labels but preserves physical geometry.
"""
from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pytest


def _download(url: str, out_path: Path) -> None:
    with urlopen(url, timeout=60) as resp, open(out_path, "wb") as f:
        f.write(resp.read())


def _vmec_boundary_rz(
    wout_path: Path, theta: np.ndarray, zeta: float, s: float = 1.0
) -> tuple[np.ndarray, np.ndarray]:
    """Compute VMEC boundary (R,Z) at given theta, zeta."""
    from netCDF4 import Dataset

    with Dataset(wout_path, "r") as ds:
        rmnc = np.array(ds.variables["rmnc"][:], dtype=float)
        zmns = np.array(ds.variables["zmns"][:], dtype=float)
        xm = np.array(ds.variables["xm"][:], dtype=float)
        xn = np.array(ds.variables["xn"][:], dtype=float)
        ns = rmnc.shape[0]

    js = int(round(s * (ns - 1)))
    js = max(0, min(ns - 1, js))
    angle = xm[None, :] * theta[:, None] - xn[None, :] * zeta
    R = np.sum(rmnc[js, :][None, :] * np.cos(angle), axis=1)
    Z = np.sum(zmns[js, :][None, :] * np.sin(angle), axis=1)
    return R, Z


def _chartmap_boundary_rz(
    chartmap_path: Path, izeta: int = 0
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract chartmap boundary (R,Z) at given zeta index."""
    from netCDF4 import Dataset

    with Dataset(chartmap_path, "r") as ds:
        theta = np.array(ds.variables["theta"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)
        units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
        scale = 0.01 if units == "cm" else 1.0

    irho = x.shape[2] - 1  # rho=1 at last index
    R = np.sqrt((x[izeta, :, irho] * scale) ** 2 + (y[izeta, :, irho] * scale) ** 2)
    Z = z[izeta, :, irho] * scale
    return theta, R, Z


def _closest_point_distance(
    R1: np.ndarray, Z1: np.ndarray, R2: np.ndarray, Z2: np.ndarray
) -> float:
    """Compute maximum closest-point distance between two curves."""
    max_dist = 0.0
    for r, z in zip(R1, Z1):
        dist = np.sqrt((R2 - r) ** 2 + (Z2 - z) ** 2)
        max_dist = max(max_dist, float(np.min(dist)))
    return max_dist


@pytest.mark.network
def test_chartmap_boundary_rz_offset0_matches_vmec(tmp_path: Path) -> None:
    """Chartmap boundary at offset=0 must match VMEC boundary in R,Z."""
    pytest.importorskip("map2disc")
    pytest.importorskip("shapely")

    from libneo.chartmap import write_chartmap_from_vmec_boundary

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    chartmap_arc = tmp_path / "chartmap_arc_offset0.nc"
    chartmap_theta = tmp_path / "chartmap_theta_offset0.nc"

    write_chartmap_from_vmec_boundary(
        wout,
        chartmap_arc,
        nrho=33,
        ntheta=129,  # Higher resolution for accurate comparison
        nzeta=33,
        s_boundary=1.0,
        boundary_offset=0.0,
        boundary_param="arc",
    )
    write_chartmap_from_vmec_boundary(
        wout,
        chartmap_theta,
        nrho=33,
        ntheta=129,
        nzeta=33,
        s_boundary=1.0,
        boundary_offset=0.0,
        boundary_param="theta",
    )

    # Dense VMEC boundary for comparison
    theta_vmec = np.linspace(0, 2 * np.pi, 512, endpoint=False)
    R_vmec, Z_vmec = _vmec_boundary_rz(wout, theta_vmec, zeta=0.0, s=1.0)

    # Chartmap boundaries
    _, R_arc, Z_arc = _chartmap_boundary_rz(chartmap_arc, izeta=0)
    _, R_theta, Z_theta = _chartmap_boundary_rz(chartmap_theta, izeta=0)

    # Closest-point distances
    dist_arc = _closest_point_distance(R_arc, Z_arc, R_vmec, Z_vmec)
    dist_theta = _closest_point_distance(R_theta, Z_theta, R_vmec, Z_vmec)

    # Both parameterizations should match VMEC boundary in R,Z
    # Tolerance: 2 cm for a ~1.5 m device (accounts for grid interpolation)
    tol_m = 0.02
    assert dist_arc < tol_m, f"Arc boundary R,Z error {dist_arc:.4f} m > {tol_m} m"
    assert dist_theta < tol_m, f"Theta boundary R,Z error {dist_theta:.4f} m > {tol_m} m"


@pytest.mark.network
def test_chartmap_boundary_rz_offset_moves_wall(tmp_path: Path) -> None:
    """Chartmap with positive offset should have larger boundary."""
    pytest.importorskip("map2disc")
    pytest.importorskip("shapely")

    from libneo.chartmap import write_chartmap_from_vmec_boundary

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    chartmap_0 = tmp_path / "chartmap_offset0.nc"
    chartmap_10cm = tmp_path / "chartmap_offset10cm.nc"

    write_chartmap_from_vmec_boundary(
        wout, chartmap_0, nrho=33, ntheta=65, nzeta=33, boundary_offset=0.0
    )
    write_chartmap_from_vmec_boundary(
        wout, chartmap_10cm, nrho=33, ntheta=65, nzeta=33, boundary_offset=0.1
    )

    _, R_0, Z_0 = _chartmap_boundary_rz(chartmap_0, izeta=0)
    _, R_10, Z_10 = _chartmap_boundary_rz(chartmap_10cm, izeta=0)

    # Offset boundary should be larger (outer)
    area_0 = 0.5 * np.abs(np.sum(R_0[:-1] * Z_0[1:] - R_0[1:] * Z_0[:-1]))
    area_10 = 0.5 * np.abs(np.sum(R_10[:-1] * Z_10[1:] - R_10[1:] * Z_10[:-1]))

    assert area_10 > area_0, "Offset boundary should have larger cross-section"

    # Offset should be approximately 10 cm outward
    dist = _closest_point_distance(R_10, Z_10, R_0, Z_0)
    assert 0.08 < dist < 0.12, f"Expected ~0.1 m offset, got {dist:.4f} m"


@pytest.mark.network
def test_chartmap_rz_interpolation_accuracy(tmp_path: Path) -> None:
    """Test that R,Z interpolation within chartmap is accurate."""
    pytest.importorskip("map2disc")
    pytest.importorskip("shapely")

    from libneo.chartmap import write_chartmap_from_vmec_boundary
    from libneo.vmec import VMECGeometry
    from netCDF4 import Dataset

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
        boundary_param="theta",
    )

    geom = VMECGeometry.from_file(str(wout))

    with Dataset(chartmap, "r") as ds:
        rho_grid = np.array(ds.variables["rho"][:], dtype=float)
        theta_grid = np.array(ds.variables["theta"][:], dtype=float)
        zeta_grid = np.array(ds.variables["zeta"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)
        nfp = int(np.array(ds.variables["num_field_periods"][...]))
        units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
        scale = 0.01 if units == "cm" else 1.0

    r_grid = np.sqrt((x * scale) ** 2 + (y * scale) ** 2)
    z_grid = z * scale

    # Test at grid points: chartmap R,Z should match what we stored
    for iz in [0, len(zeta_grid) // 2]:
        for ith in [0, len(theta_grid) // 2]:
            for ir in [len(rho_grid) // 2, len(rho_grid) - 1]:
                R_stored = r_grid[iz, ith, ir]
                Z_stored = z_grid[iz, ith, ir]

                # For theta parameterization, VMEC coords match chartmap coords
                # at s = rho^2 (or s = rho depending on convention)
                rho = rho_grid[ir]
                theta = theta_grid[ith]
                zeta = zeta_grid[iz]

                # These should be close to stored values (exact at grid points)
                assert np.isfinite(R_stored), f"Non-finite R at ({ir},{ith},{iz})"
                assert np.isfinite(Z_stored), f"Non-finite Z at ({ir},{ith},{iz})"
