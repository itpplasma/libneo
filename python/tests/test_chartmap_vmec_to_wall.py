"""Tests for VMEC chartmap with extension to wall boundaries."""
from __future__ import annotations

from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pytest


def _download(url: str, out_path: Path) -> None:
    with urlopen(url, timeout=60) as resp, open(out_path, "wb") as f:
        f.write(resp.read())


def _figure_path(name: str) -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    target = repo_root / "build" / "test" / "python"
    target.mkdir(parents=True, exist_ok=True)
    return target / name


def _plot_vmec_to_wall_grid(
    chartmap_path: Path,
    wout_path: Path,
    out_path: Path,
    ir_lcfs: int,
    wall_rz: tuple[np.ndarray, np.ndarray] | None = None,
    izeta: int = 0,
) -> None:
    """Plot R,Z grid showing VMEC interior (blue) and wall extension (orange)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from netCDF4 import Dataset

    from libneo.vmec import VMECGeometry

    with Dataset(chartmap_path, "r") as ds:
        rho_grid = np.array(ds.variables["rho"][:], dtype=float)
        theta_grid = np.array(ds.variables["theta"][:], dtype=float)
        zeta_grid = np.array(ds.variables["zeta"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)
        units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
        scale = 0.01 if units == "cm" else 1.0

    R = np.sqrt((x * scale) ** 2 + (y * scale) ** 2)
    Z = z * scale

    geom = VMECGeometry.from_file(str(wout_path))

    fig, ax = plt.subplots(figsize=(8.0, 8.0))

    zeta = float(zeta_grid[izeta])
    nrho = len(rho_grid)
    ntheta = len(theta_grid)

    for ir in range(nrho):
        color = "C0" if ir <= ir_lcfs else "C1"
        lw = 1.5 if ir == ir_lcfs or ir == nrho - 1 else 0.6
        r_ring = R[izeta, :, ir]
        z_ring = Z[izeta, :, ir]
        ax.plot(
            np.r_[r_ring, r_ring[0]],
            np.r_[z_ring, z_ring[0]],
            color=color,
            linewidth=lw,
            alpha=0.9 if ir == ir_lcfs or ir == nrho - 1 else 0.7,
        )

    for ith in range(0, ntheta, max(1, ntheta // 24)):
        ax.plot(
            R[izeta, ith, : ir_lcfs + 1],
            Z[izeta, ith, : ir_lcfs + 1],
            color="C0",
            linewidth=0.4,
            alpha=0.5,
        )
        ax.plot(
            R[izeta, ith, ir_lcfs:],
            Z[izeta, ith, ir_lcfs:],
            color="C1",
            linewidth=0.4,
            alpha=0.5,
        )

    theta_dense = np.linspace(0.0, 2.0 * np.pi, 256, endpoint=False)
    R_lcfs, Z_lcfs, _ = geom.coords_s(1.0, theta_dense, zeta, use_asym=True)
    ax.plot(R_lcfs, Z_lcfs, "k--", linewidth=1.0, label="VMEC LCFS")

    if wall_rz is not None:
        R_wall, Z_wall = wall_rz
        ax.plot(
            np.r_[R_wall, R_wall[0]],
            np.r_[Z_wall, Z_wall[0]],
            "r-",
            linewidth=1.5,
            label="Wall boundary",
        )

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title(f"VMEC + wall extension (zeta={zeta:.3f})")

    legend_elements = [
        Line2D([0], [0], color="C0", linewidth=2, label="VMEC region"),
        Line2D([0], [0], color="C1", linewidth=2, label="Extension to wall"),
        Line2D([0], [0], color="k", linewidth=1, linestyle="--", label="VMEC LCFS"),
    ]
    if wall_rz is not None:
        legend_elements.append(
            Line2D([0], [0], color="r", linewidth=1.5, label="Wall boundary")
        )
    ax.legend(handles=legend_elements, loc="upper right")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def _create_synthetic_wall(
    geom, zeta: float, offset: float, ntheta: int = 512
) -> tuple[np.ndarray, np.ndarray]:
    """Create a synthetic wall boundary using shapely buffer (handles concave regions)."""
    from shapely.geometry import Polygon

    theta = np.linspace(0.0, 2.0 * np.pi, 4096, endpoint=False)
    R_lcfs, Z_lcfs, _ = geom.coords_s(1.0, theta, zeta, use_asym=True)

    poly = Polygon(np.column_stack([R_lcfs, Z_lcfs]))
    if not poly.is_valid:
        poly = poly.buffer(0.0)

    off = poly.buffer(offset, join_style=1)
    if off.geom_type != "Polygon":
        off = max(list(off.geoms), key=lambda g: g.area)

    coords = np.asarray(off.exterior.coords, dtype=float)
    coords_open = coords[:-1, :]

    i0 = int(np.argmax(coords_open[:, 0]))
    coords_open = np.vstack([coords_open[i0:, :], coords_open[:i0, :]])
    if coords_open.shape[0] >= 2 and coords_open[1, 1] - coords_open[0, 1] < 0.0:
        coords_open = coords_open[::-1, :]

    seg = np.sqrt(np.sum(np.diff(coords_open, axis=0) ** 2, axis=1))
    s_coords = np.concatenate(([0.0], np.cumsum(seg)))
    length = s_coords[-1]

    s_query = np.linspace(0.0, length, ntheta, endpoint=False)
    R_wall = np.interp(s_query, s_coords, coords_open[:, 0])
    Z_wall = np.interp(s_query, s_coords, coords_open[:, 1])

    return R_wall, Z_wall


@pytest.mark.network
def test_chartmap_vmec_to_wall_synthetic(tmp_path: Path) -> None:
    """Test VMEC to wall mode with synthetic wall boundary."""
    from netCDF4 import Dataset

    from libneo.chartmap import write_chartmap_from_vmec_to_wall
    from libneo.vmec import VMECGeometry

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    geom = VMECGeometry.from_file(str(wout))

    with Dataset(wout, "r") as ds:
        nfp = int(np.array(ds.variables["nfp"][...]))

    nrho = 33
    ntheta = 65
    nzeta = 5
    rho_lcfs = 0.8
    wall_offset = 0.15

    period = 2.0 * np.pi / nfp
    wall_zeta = np.linspace(0.0, period, nzeta, endpoint=False, dtype=float)
    wall_rz = []
    for zeta in wall_zeta:
        R_wall, Z_wall = _create_synthetic_wall(geom, float(zeta), wall_offset)
        wall_rz.append((R_wall, Z_wall))

    chartmap = tmp_path / "chartmap_vmec_to_wall.nc"
    write_chartmap_from_vmec_to_wall(
        wout,
        chartmap,
        wall_rz,
        wall_zeta,
        nrho=nrho,
        ntheta=ntheta,
        rho_lcfs=rho_lcfs,
        num_field_periods=nfp,
    )

    with Dataset(chartmap, "r") as ds:
        rho_grid = np.array(ds.variables["rho"][:], dtype=float)
        theta_grid = np.array(ds.variables["theta"][:], dtype=float)
        zeta_grid = np.array(ds.variables["zeta"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)
        rho_conv = str(ds.getncattr("rho_convention"))
        units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
        scale = 0.01 if units == "cm" else 1.0

    assert rho_conv == "vmec_to_wall"
    assert len(zeta_grid) == nzeta
    assert np.allclose(zeta_grid, wall_zeta)

    R = np.sqrt((x * scale) ** 2 + (y * scale) ** 2)
    Z = z * scale

    ir_lcfs = int(np.searchsorted(rho_grid, rho_lcfs, side="right") - 1)
    rho_lcfs_actual = float(rho_grid[ir_lcfs])

    for iz in [0, len(zeta_grid) // 2]:
        zeta = float(zeta_grid[iz])
        theta = theta_grid

        for ir in range(ir_lcfs + 1):
            rho_val = rho_grid[ir]
            s_val = (float(rho_val) / rho_lcfs_actual) ** 2
            R_vmec, Z_vmec, _ = geom.coords_s(s_val, theta, zeta, use_asym=True)
            assert np.allclose(R_vmec, R[iz, :, ir], rtol=0.0, atol=1.0e-6), (
                f"VMEC region mismatch at iz={iz}, ir={ir}"
            )
            assert np.allclose(Z_vmec, Z[iz, :, ir], rtol=0.0, atol=1.0e-6), (
                f"VMEC region mismatch at iz={iz}, ir={ir}"
            )

    R_outer = R[0, :, -1]
    Z_outer = Z[0, :, -1]
    R_wall_expected, Z_wall_expected = wall_rz[0]
    d2 = (R_outer[:, None] - R_wall_expected[None, :]) ** 2 + (
        Z_outer[:, None] - Z_wall_expected[None, :]
    ) ** 2
    min_dist = np.sqrt(np.min(d2, axis=1))
    max_dist = float(np.max(min_dist))
    assert max_dist < 5.0e-3, f"outer boundary too far from wall (max dist {max_dist:.2e} m)"


@pytest.mark.network
def test_chartmap_vmec_to_wall_poloidal_normal_matching(tmp_path: Path) -> None:
    """Wall matching along LCFS poloidal-plane normal should preserve offsets."""
    from netCDF4 import Dataset

    from libneo.chartmap import write_chartmap_from_vmec_to_wall
    from libneo.vmec import VMECGeometry

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    geom = VMECGeometry.from_file(str(wout))

    with Dataset(wout, "r") as ds:
        nfp = int(np.array(ds.variables["nfp"][...]))

    nrho = 33
    ntheta = 65
    nzeta = 5
    rho_lcfs = 0.8
    wall_offset = 0.15

    period = 2.0 * np.pi / nfp
    wall_zeta = np.linspace(0.0, period, nzeta, endpoint=False, dtype=float)

    wall_rz = []
    theta_wall = np.linspace(0.0, 2.0 * np.pi, 1024, endpoint=False, dtype=float)
    dtheta_wall = float(theta_wall[1] - theta_wall[0])

    for zeta in wall_zeta:
        R_lcfs, Z_lcfs, dR_ds, dZ_ds = geom.coords_s_with_deriv(
            1.0, theta_wall, float(zeta), use_asym=True
        )

        dR_dth = (np.roll(R_lcfs, -1) - np.roll(R_lcfs, 1)) / (2.0 * dtheta_wall)
        dZ_dth = (np.roll(Z_lcfs, -1) - np.roll(Z_lcfs, 1)) / (2.0 * dtheta_wall)

        nR = dZ_dth
        nZ = -dR_dth

        nmag = np.sqrt(nR**2 + nZ**2)
        nmag = np.where(nmag == 0.0, 1.0, nmag)
        nR = nR / nmag
        nZ = nZ / nmag

        outward_dot = nR * dR_ds + nZ * dZ_ds
        flip = outward_dot < 0.0
        nR = np.where(flip, -nR, nR)
        nZ = np.where(flip, -nZ, nZ)

        R_wall = R_lcfs + wall_offset * nR
        Z_wall = Z_lcfs + wall_offset * nZ
        wall_rz.append((R_wall, Z_wall))

    chartmap = tmp_path / "chartmap_vmec_to_wall_normal.nc"
    write_chartmap_from_vmec_to_wall(
        wout,
        chartmap,
        wall_rz,
        wall_zeta,
        nrho=nrho,
        ntheta=ntheta,
        rho_lcfs=rho_lcfs,
        num_field_periods=nfp,
        wall_match="poloidal_normal",
    )

    with Dataset(chartmap, "r") as ds:
        rho_grid = np.array(ds.variables["rho"][:], dtype=float)
        theta_grid = np.array(ds.variables["theta"][:], dtype=float)
        zeta_grid = np.array(ds.variables["zeta"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)
        units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
        scale = 0.01 if units == "cm" else 1.0

    R = np.sqrt((x * scale) ** 2 + (y * scale) ** 2)
    Z = z * scale

    ir_lcfs = int(np.searchsorted(rho_grid, rho_lcfs, side="right") - 1)
    rho_lcfs_actual = float(rho_grid[ir_lcfs])

    iz = 0
    zeta0 = float(zeta_grid[iz])
    R_lcfs_g, Z_lcfs_g, dR_ds_g, dZ_ds_g = geom.coords_s_with_deriv(
        1.0, theta_grid, zeta0, use_asym=True
    )

    dtheta = float(theta_grid[1] - theta_grid[0])
    dR_dth = (np.roll(R_lcfs_g, -1) - np.roll(R_lcfs_g, 1)) / (2.0 * dtheta)
    dZ_dth = (np.roll(Z_lcfs_g, -1) - np.roll(Z_lcfs_g, 1)) / (2.0 * dtheta)

    nR = dZ_dth
    nZ = -dR_dth
    nmag = np.sqrt(nR**2 + nZ**2)
    nmag = np.where(nmag == 0.0, 1.0, nmag)
    nR = nR / nmag
    nZ = nZ / nmag
    outward_dot = nR * dR_ds_g + nZ * dZ_ds_g
    flip = outward_dot < 0.0
    nR = np.where(flip, -nR, nR)
    nZ = np.where(flip, -nZ, nZ)

    R_expected = R_lcfs_g + wall_offset * nR
    Z_expected = Z_lcfs_g + wall_offset * nZ

    assert np.allclose(R[iz, :, -1], R_expected, rtol=0.0, atol=5.0e-4)
    assert np.allclose(Z[iz, :, -1], Z_expected, rtol=0.0, atol=5.0e-4)

    _plot_vmec_to_wall_grid(
        chartmap,
        wout,
        _figure_path("chartmap_vmec_to_wall_ncsx.png"),
        ir_lcfs=ir_lcfs,
        wall_rz=wall_rz[0],
        izeta=0,
    )


@pytest.mark.network
def test_chartmap_vmec_to_wall_c1_continuity(tmp_path: Path) -> None:
    """Test that derivatives are continuous at LCFS (C1 continuity)."""
    from netCDF4 import Dataset

    from libneo.chartmap import write_chartmap_from_vmec_to_wall
    from libneo.vmec import VMECGeometry

    wout_url = (
        "https://princetonuniversity.github.io/STELLOPT/examples/wout_ncsx_c09r00_fixed.nc"
    )
    wout = tmp_path / "wout_ncsx.nc"
    _download(wout_url, wout)

    geom = VMECGeometry.from_file(str(wout))

    with Dataset(wout, "r") as ds:
        nfp = int(np.array(ds.variables["nfp"][...]))

    nrho = 65
    ntheta = 65
    nzeta = 5
    rho_lcfs = 0.75
    wall_offset = 0.12

    period = 2.0 * np.pi / nfp
    wall_zeta = np.linspace(0.0, period, nzeta, endpoint=False, dtype=float)
    wall_rz = []
    for zeta in wall_zeta:
        R_wall, Z_wall = _create_synthetic_wall(geom, float(zeta), wall_offset)
        wall_rz.append((R_wall, Z_wall))

    chartmap = tmp_path / "chartmap_continuity.nc"
    write_chartmap_from_vmec_to_wall(
        wout,
        chartmap,
        wall_rz,
        wall_zeta,
        nrho=nrho,
        ntheta=ntheta,
        rho_lcfs=rho_lcfs,
        num_field_periods=nfp,
    )

    with Dataset(chartmap, "r") as ds:
        rho_grid = np.array(ds.variables["rho"][:], dtype=float)
        x = np.array(ds.variables["x"][:], dtype=float)
        y = np.array(ds.variables["y"][:], dtype=float)
        z = np.array(ds.variables["z"][:], dtype=float)
        units = str(getattr(ds.variables["x"], "units", "cm")).strip().lower()
        scale = 0.01 if units == "cm" else 1.0

    R = np.sqrt((x * scale) ** 2 + (y * scale) ** 2)
    Z = z * scale

    ir_lcfs = int(np.searchsorted(rho_grid, rho_lcfs, side="right") - 1)

    for iz in range(Z.shape[0]):
        for ith in range(0, Z.shape[1], 8):
            R_line = R[iz, ith, :]
            Z_line = Z[iz, ith, :]

            dR_drho = np.gradient(R_line, rho_grid)
            dZ_drho = np.gradient(Z_line, rho_grid)

            dR_jump = abs(dR_drho[ir_lcfs + 1] - dR_drho[ir_lcfs])
            dZ_jump = abs(dZ_drho[ir_lcfs + 1] - dZ_drho[ir_lcfs])

            dR_scale = max(abs(dR_drho[ir_lcfs]), abs(dR_drho[ir_lcfs + 1]), 0.01)
            dZ_scale = max(abs(dZ_drho[ir_lcfs]), abs(dZ_drho[ir_lcfs + 1]), 0.01)

            assert dR_jump / dR_scale < 2.0, (
                f"dR/drho discontinuity at iz={iz}, ith={ith}: jump={dR_jump:.4f}, "
                f"scale={dR_scale:.4f}, ratio={dR_jump/dR_scale:.2f}"
            )
            assert dZ_jump / dZ_scale < 2.0, (
                f"dZ/drho discontinuity at iz={iz}, ith={ith}: jump={dZ_jump:.4f}, "
                f"scale={dZ_scale:.4f}, ratio={dZ_jump/dZ_scale:.2f}"
            )
